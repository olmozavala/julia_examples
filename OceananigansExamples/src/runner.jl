const RUNS_DIR = abspath(joinpath(@__DIR__, "..", "runs"))
const NETCDF_FILENAME = "simulation.nc"
const RUNNING_MARKER_FILENAME = "running.lock"
const COMPLETED_MARKER_FILENAME = "completed.ok"
const FAILED_MARKER_FILENAME = "failed.txt"

Base.@kwdef mutable struct SimulationSession
    id::String = "idle"
    example_id::String = "barotropic_jet"
    run_name::String = "idle"
    run_dir::String = RUNS_DIR
    netcdf_path::String = ""
    requested_backend::Symbol = :gpu
    effective_backend::Symbol = :cpu
    backend_label::String = "CPU"
    backend_note::String = "Idle"
    lock::ReentrantLock = ReentrantLock()
    running::Bool = false
    completed::Bool = false
    failed::Bool = false
    status::String = "Idle"
    error_message::String = ""
    stop_time::Float64 = 0.0
    output_interval::Int = 1
    started_at::DateTime = now()
    completed_at::Union{Nothing, DateTime} = nothing
    task::Union{Task, Nothing} = nothing
end

const RUN_SESSIONS = SimulationSession[]
const RUN_SESSIONS_LOCK = ReentrantLock()

"""
    centered_along_x(raw::AbstractArray{<:Real, 3}, Nx::Int)

Centers the `raw` data along the x-axis by averaging adjacent cells in the x-direction.
"""
function centered_along_x(raw::AbstractArray{<:Real, 3}, Nx::Int)
    size(raw, 1) == Nx + 1 || return Array(raw)
    return @views 0.5 .* (raw[1:Nx, :, :] .+ raw[2:Nx+1, :, :])
end

"""
    centered_along_y(raw::AbstractArray{<:Real, 3}, Ny::Int)

Centers the `raw` data along the y-axis by averaging adjacent cells in the y-direction.
"""
function centered_along_y(raw::AbstractArray{<:Real, 3}, Ny::Int)
    size(raw, 2) == Ny + 1 || return Array(raw)
    return @views 0.5 .* (raw[:, 1:Ny, :] .+ raw[:, 2:Ny+1, :])
end

"""
    depth_representative(field::Array{Float64, 3}; mode::Symbol = :mean)

Calculates the depth representative value of `field` either by taking the mean across the depth (if `mode = :mean`) or returning the surface slice (if `mode = :surface`).
"""
function depth_representative(field::Array{Float64, 3}; mode::Symbol = :mean)
    mode == :surface && return Array(field[:, :, end])
    return dropdims(mean(field; dims = 3), dims = 3)
end

"""
    finite_difference_vorticity(u::Matrix{Float64}, v::Matrix{Float64}, dx::Float64, dy::Float64)

Computes the relative vorticity from velocity fields `u` and `v` using a finite difference scheme over grid points separated by `dx` and `dy`.
"""
function finite_difference_vorticity(u::Matrix{Float64}, v::Matrix{Float64}, dx::Float64, dy::Float64)
    Nx, Ny = size(u)
    ζ = zeros(Float64, Nx, Ny)
    @inbounds for j in 1:Ny
        jp = j == Ny ? 1 : j + 1
        jm = j == 1 ? Ny : j - 1
        for i in 1:Nx
            ip = i == Nx ? 1 : i + 1
            im = i == 1 ? Nx : i - 1
            dvdx = (v[ip, j] - v[im, j]) / (2dx)
            dudy = (u[i, jp] - u[i, jm]) / (2dy)
            ζ[i, j] = dvdx - dudy
        end
    end
    return ζ
end

"""
    field_outputs(model)

Extracts and returns the velocity and tracer fields from the simulator `model` as a dictionary.
"""
function field_outputs(model)
    outputs = Dict{String, Any}()
    for name in propertynames(model.velocities)
        outputs[string(name)] = getproperty(model.velocities, name)
    end
    for name in propertynames(model.tracers)
        outputs[string(name)] = getproperty(model.tracers, name)
    end
    return outputs
end

"""
    depth_mean_speed_data(model, settings::RunSettings)

Computes the depth-mean speed data from the `model` velocities after centering the velocities onto the correct grid nodes.
"""
function depth_mean_speed_data(model, settings::RunSettings)
    d = settings.domain
    raw_u = Array(interior(model.velocities.u))
    raw_v = Array(interior(model.velocities.v))
    u = centered_along_x(raw_u, d.Nx)
    v = centered_along_y(raw_v, d.Ny)
    u_bar = depth_representative(u)
    v_bar = depth_representative(v)
    return sqrt.(u_bar .^ 2 .+ v_bar .^ 2)
end

"""
    depth_mean_relative_vorticity_data(model, settings::RunSettings)

Computes the depth-averaged relative vorticity data by extracting the velocities and using finite difference operators.
"""
function depth_mean_relative_vorticity_data(model, settings::RunSettings)
    d = settings.domain
    raw_u = Array(interior(model.velocities.u))
    raw_v = Array(interior(model.velocities.v))
    u = centered_along_x(raw_u, d.Nx)
    v = centered_along_y(raw_v, d.Ny)
    u_bar = depth_representative(u)
    v_bar = depth_representative(v)
    dx = d.Lx / d.Nx
    dy = d.Ly / d.Ny
    return finite_difference_vorticity(u_bar, v_bar, dx, dy)
end

"""
    kinetic_energy_data(model, settings::RunSettings)

Calculates the domain-averaged kinetic energy.
"""
function kinetic_energy_data(model, settings::RunSettings)
    d = settings.domain
    raw_u = Array(interior(model.velocities.u))
    raw_v = Array(interior(model.velocities.v))
    u = centered_along_x(raw_u, d.Nx)
    v = centered_along_y(raw_v, d.Ny)
    return 0.5 * mean(u .^ 2 .+ v .^ 2)
end

"""
    enstrophy_data(model, settings::RunSettings)

Calculates the grid's enstrophy from the mean relative vorticity squared.
"""
enstrophy_data(model, settings::RunSettings) = mean(depth_mean_relative_vorticity_data(model, settings) .^ 2)
"""
    max_speed_data(model, settings::RunSettings)

Calculates the maximum local speed in the domain.
"""
max_speed_data(model, settings::RunSettings) = maximum(depth_mean_speed_data(model, settings))

"""
    derived_outputs(model, settings::RunSettings)

Packages diagnostic and derived values (e.g. vorticity, speed, kinetic energy) as dictionary closures for exporting to file.
"""
function derived_outputs(model, settings::RunSettings)
    outputs = Dict{String, Any}(
        "depth_mean_relative_vorticity" => (m -> depth_mean_relative_vorticity_data(m, settings)),
        "depth_mean_speed" => (m -> depth_mean_speed_data(m, settings)),
        "kinetic_energy" => (m -> kinetic_energy_data(m, settings)),
        "enstrophy" => (m -> enstrophy_data(m, settings)),
        "max_speed" => (m -> max_speed_data(m, settings)),
    )

    if :b in propertynames(model.tracers)
        outputs["surface_buoyancy"] = function (m)
            raw_b = Array(interior(m.tracers.b))
            return depth_representative(raw_b; mode = :surface)
        end
    end

    return outputs
end

"""
    output_dimensions(model)

Defines the grid dimensionality for each specific output field generated by `field_outputs` and `derived_outputs`.
"""
function output_dimensions(model)
    dims = Dict{String, Tuple}()
    dims["depth_mean_relative_vorticity"] = ("x_caa", "y_aca")
    dims["depth_mean_speed"] = ("x_caa", "y_aca")
    dims["kinetic_energy"] = ()
    dims["enstrophy"] = ()
    dims["max_speed"] = ()

    if :b in propertynames(model.tracers)
        dims["surface_buoyancy"] = ("x_caa", "y_aca")
    end

    return dims
end

"""
    output_attributes()

Constructs the output attributes (like units or long names) used in the NetCDF output variables.
"""
function output_attributes()
    return Dict(
        "depth_mean_relative_vorticity" => Dict("long_name" => "Depth-mean relative vorticity"),
        "depth_mean_speed" => Dict("long_name" => "Depth-mean speed", "units" => "m s-1"),
        "surface_buoyancy" => Dict("long_name" => "Surface buoyancy"),
        "kinetic_energy" => Dict("long_name" => "Domain-mean kinetic energy"),
        "enstrophy" => Dict("long_name" => "Domain-mean enstrophy"),
        "max_speed" => Dict("long_name" => "Depth-mean maximum speed", "units" => "m s-1"),
    )
end

"""
    unique_run_stamp()

Generates a timestamped string to create a uniquely identifiable run identifier.
"""
function unique_run_stamp()
    return Dates.format(now(), dateformat"yyyymmdd_HHMMSS_s")
end

"""
    create_run_directory(example_id::String)

Sets up and allocates an output directory for a specific `example_id` run, returning its path.
"""
function create_run_directory(example_id::String)
    mkpath(RUNS_DIR)
    base = "$(example_id)_$(unique_run_stamp())"
    run_dir = joinpath(RUNS_DIR, base)
    suffix = 1
    while ispath(run_dir)
        suffix += 1
        run_dir = joinpath(RUNS_DIR, "$(base)_$(suffix)")
    end
    mkpath(run_dir)
    return run_dir
end

"""
    run_marker_path(run_dir::String, filename::String)

Computes the file path to a given status marker inside a run output directory.
"""
run_marker_path(run_dir::String, filename::String) = joinpath(run_dir, filename)

"""
    write_run_marker(run_dir::String, filename::String, message::String)

Writes the string `message` into a specified `filename` inside `run_dir` to act as a lock or status tracking file.
"""
function write_run_marker(run_dir::String, filename::String, message::String)
    open(run_marker_path(run_dir, filename), "w") do io
        write(io, message)
        endswith(message, "\n") || write(io, "\n")
    end
    return nothing
end

"""
    remove_run_marker(run_dir::String, filename::String)

Purges the specified marker file from `run_dir` if it exists.
"""
function remove_run_marker(run_dir::String, filename::String)
    path = run_marker_path(run_dir, filename)
    isfile(path) && rm(path; force = true)
    return nothing
end

"""
    register_session!(session::SimulationSession)

Threads-safely appends a new `SimulationSession` to the global `RUN_SESSIONS` list.
"""
function register_session!(session::SimulationSession)
    lock(RUN_SESSIONS_LOCK) do
        push!(RUN_SESSIONS, session)
    end
    return session
end

"""
    list_sessions()

Returns a safely sorted copy of the current `RUN_SESSIONS`, prioritized by started time.
"""
function list_sessions()
    lock(RUN_SESSIONS_LOCK) do
        return sort(copy(RUN_SESSIONS); by = session -> session.started_at, rev = true)
    end
end

"""
    update_status!(session::SimulationSession, simulation)

Updates the user-facing status string for a running `simulation` inside a `session` to report time and iteration values.
"""
function update_status!(session::SimulationSession, simulation)
    lock(session.lock) do
        session.status = @sprintf(
            "Running %s on %s | t = %.2f h | iter = %d",
            replace(session.example_id, '_' => ' '),
            session.backend_label,
            simulation.model.clock.time / 3600,
            simulation.model.clock.iteration,
        )
    end
    return nothing
end

"""
    finalize_session!(session::SimulationSession; status="Completed", failed=false, error_message="")

Flags a simulation session as completed or failed and writes the corresponding completion or fail marker files.
"""
function finalize_session!(session::SimulationSession; status = "Completed", failed = false, error_message = "")
    completed_at = now()
    remove_run_marker(session.run_dir, RUNNING_MARKER_FILENAME)
    remove_run_marker(session.run_dir, failed ? COMPLETED_MARKER_FILENAME : FAILED_MARKER_FILENAME)

    if failed
        write_run_marker(
            session.run_dir,
            FAILED_MARKER_FILENAME,
            "completed_at=$(completed_at)\nstatus=$(status)\nerror=$(error_message)",
        )
    else
        write_run_marker(
            session.run_dir,
            COMPLETED_MARKER_FILENAME,
            "completed_at=$(completed_at)\nstatus=$(status)\nnetcdf=$(session.netcdf_path)",
        )
    end

    lock(session.lock) do
        session.running = false
        session.completed = !failed
        session.failed = failed
        session.status = status
        session.error_message = error_message
        session.completed_at = completed_at
    end
    return nothing
end

"""
    run_global_attributes(session::SimulationSession, settings::RunSettings, spec::ExampleSpec)

Creates a descriptive dictionary of global properties (like domain, dt, simulation setup) saved in the NetCDF archive.
"""
function run_global_attributes(session::SimulationSession, settings::RunSettings, spec::ExampleSpec)
    d = settings.domain
    return Dict(
        "run_name" => session.run_name,
        "example_id" => spec.id,
        "example_label" => spec.label,
        "description" => spec.description,
        "requested_backend" => String(session.requested_backend),
        "effective_backend" => String(session.effective_backend),
        "backend_label" => session.backend_label,
        "backend_note" => session.backend_note,
        "started_at" => string(session.started_at),
        "dt" => settings.dt,
        "stop_time" => settings.stop_time,
        "output_interval" => settings.output_interval,
        "Lx" => d.Lx,
        "Ly" => d.Ly,
        "Lz" => d.Lz,
        "Nx" => d.Nx,
        "Ny" => d.Ny,
        "Nz" => d.Nz,
    )
end

"""
    run_loop!(session::SimulationSession, settings::RunSettings, backend)

Initializes the Oceananigans internal model, sets up NetCDF writers with scheduled output cycles, and runs the simulation loop.
"""
function run_loop!(session::SimulationSession, settings::RunSettings, backend)
    spec = EXAMPLE_SPECS[settings.example_id]
    model = spec.builder(settings, backend.arch)
    simulation = Simulation(model; Δt = settings.dt, stop_time = settings.stop_time)

    outputs = merge(field_outputs(model), derived_outputs(model, settings))
    simulation.output_writers[:netcdf] = NetCDFWriter(
        model,
        outputs;
        filename = session.netcdf_path,
        schedule = IterationInterval(settings.output_interval),
        with_halos = false,
        include_grid_metrics = false,
        dimensions = output_dimensions(model),
        output_attributes = output_attributes(),
        global_attributes = run_global_attributes(session, settings, spec),
    )

    add_callback!(simulation, sim -> update_status!(session, sim), IterationInterval(settings.output_interval); name = :progress)
    update_status!(session, simulation)
    run!(simulation)

    return finalize_session!(
        session;
        status = @sprintf("Completed %s | %s", session.run_name, basename(session.netcdf_path)),
    )
end

"""
    start_run!(settings::RunSettings)

Registers a new background `SimulationSession` via `register_session!`, creates a trackable run directory, and fires off `run_loop!` in an asynchronous worker thread.
"""
function start_run!(settings::RunSettings)
    backend = resolve_architecture(settings.preferred_backend)
    run_dir = create_run_directory(settings.example_id)
    run_name = basename(run_dir)
    netcdf_path = joinpath(run_dir, NETCDF_FILENAME)
    write_run_marker(
        run_dir,
        RUNNING_MARKER_FILENAME,
        "started_at=$(now())\nexample_id=$(settings.example_id)\nbackend=$(settings.preferred_backend)\nnetcdf=$(netcdf_path)",
    )

    session = register_session!(SimulationSession(
        id = run_name,
        example_id = settings.example_id,
        run_name = run_name,
        run_dir = run_dir,
        netcdf_path = netcdf_path,
        requested_backend = settings.preferred_backend,
        effective_backend = backend.effective,
        backend_label = backend.label,
        backend_note = backend.note,
        running = true,
        stop_time = settings.stop_time,
        output_interval = settings.output_interval,
        started_at = now(),
        status = "Building Oceananigans model for $(run_name)…",
    ))

    session.task = Threads.@spawn begin
        try
            run_loop!(session, settings, backend)
        catch err
            finalize_session!(
                session;
                status = "Run failed: $(session.run_name)",
                failed = true,
                error_message = sprint(showerror, err),
            )
        end
    end

    return session
end
