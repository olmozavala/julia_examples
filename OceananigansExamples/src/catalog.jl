const DEFAULT_FIELD_OPTIONS = (
    (label = "Relative vorticity", value = "relative_vorticity"),
    (label = "Speed", value = "speed"),
    (label = "Buoyancy", value = "buoyancy"),
)

const DEFAULT_BACKEND_OPTIONS = (
    (label = "GPU (default)", value = "gpu"),
    (label = "CPU", value = "cpu"),
)

Base.@kwdef struct DomainSettings
    Lx::Float64 = 1_200e3
    Ly::Float64 = 900e3
    Lz::Float64 = 400.0
    Nx::Int = 128
    Ny::Int = 96
    Nz::Int = 12
end

Base.@kwdef struct RunSettings
    example_id::String
    dt::Float64
    stop_time::Float64
    output_interval::Int
    domain::DomainSettings
    display_variable::Symbol
    preferred_backend::Symbol = :gpu
end

struct ExampleSpec
    id::String
    label::String
    description::String
    default_dt::Float64
    default_stop_time::Float64
    default_output_interval::Int
    default_domain::DomainSettings
    default_display_variable::Symbol
    builder::Function
end

"""
    supports_buoyancy(spec::ExampleSpec)

Checks if the target `ExampleSpec` supports and displays buoyancy.
"""
supports_buoyancy(spec::ExampleSpec) = spec.default_display_variable == :buoyancy

"""
    sech(x)

Computes the hyperbolic secant of `x`.
"""
sech(x) = inv(cosh(x))

"""
    gaussian_vortex_streamfunction(x, y, x0, y0, radius, amplitude)

Returns the streamfunction value for a localized gaussian vortex centered at `(x0, y0)`.
"""
gaussian_vortex_streamfunction(x, y, x0, y0, radius, amplitude) =
    amplitude * exp(-((x - x0)^2 + (y - y0)^2) / radius^2)

"""
    gaussian_vortex_velocity(x, y, x0, y0, radius, amplitude)

Returns the velocity components (u, v) and streamfunction derived for a 2D gaussian vortex centered at `(x0, y0)`.
"""
function gaussian_vortex_velocity(x, y, x0, y0, radius, amplitude)
    ψ = gaussian_vortex_streamfunction(x, y, x0, y0, radius, amplitude)
    u = 2 * amplitude * (y - y0) / radius^2 * exp(-((x - x0)^2 + (y - y0)^2) / radius^2)
    v = -2 * amplitude * (x - x0) / radius^2 * exp(-((x - x0)^2 + (y - y0)^2) / radius^2)
    return u, v, ψ
end

"""
    resolve_architecture(preferred_backend::Symbol)

Checks if the desired `preferred_backend` translates to a valid CUDA architecture or rolls back to `CPU()`.
"""
function resolve_architecture(preferred_backend::Symbol)
    if preferred_backend == :cpu
        return (arch = CPU(), effective = :cpu, label = "CPU", note = "Running on CPU by request.")
    end

    if CUDA.functional() && CUDA.has_cuda_gpu()
        return (arch = GPU(), effective = :gpu, label = "GPU", note = "Running on GPU by default.")
    end

    return (
        arch = CPU(),
        effective = :cpu,
        label = "CPU fallback",
        note = "CUDA is unavailable, so the run fell back to CPU.",
    )
end

"""
    periodic_model(settings::RunSettings, arch; ν, κ=ν, tracers=nothing, buoyancy=nothing)

Instantiates a Rectilinear `NonhydrostaticModel` domain parameterized with periodic boundary conditions along x/y and bounds across z.
"""
function periodic_model(settings::RunSettings, arch; ν, κ = ν, tracers = nothing, buoyancy = nothing)
    d = settings.domain
    grid = RectilinearGrid(
        arch,
        size = (d.Nx, d.Ny, d.Nz),
        x = (0, d.Lx),
        y = (0, d.Ly),
        z = (-d.Lz, 0),
        topology = (Periodic, Periodic, Bounded),
    )
    kwargs = (; grid, closure = ScalarDiffusivity(ν = ν, κ = κ))
    if tracers === nothing
        return NonhydrostaticModel(; kwargs...)
    end
    return NonhydrostaticModel(; kwargs..., tracers, buoyancy)
end

"""
    barotropic_jet_builder(settings::RunSettings, arch)

Provides an initial condition specification representing a barotropic jet meandering slightly across the given domain setup.
"""
function barotropic_jet_builder(settings::RunSettings, arch)
    d = settings.domain
    width = 0.10 * d.Ly
    perturbation_width = 0.08 * d.Ly
    U0 = 0.9
    model = periodic_model(settings, arch; ν = 20.0)
    u0(x, y, z) = U0 * tanh((y - 0.5d.Ly) / width)
    v0(x, y, z) = 0.08U0 * sin(4π * x / d.Lx) * exp(-((y - 0.5d.Ly)^2) / perturbation_width^2)
    set!(model, u = u0, v = v0)
    return model
end

"""
    baroclinic_jet_builder(settings::RunSettings, arch)

Builds a baroclinic structure jet with prescribed buoyancy gradients that is surface-intensified.
"""
function baroclinic_jet_builder(settings::RunSettings, arch)
    d = settings.domain
    width = 0.09 * d.Ly
    front_width = 0.12 * d.Ly
    U0 = 0.65
    N² = 4e-5
    B0 = 1.8e-3
    model = periodic_model(settings, arch; ν = 10.0, κ = 5.0, tracers = :b, buoyancy = BuoyancyTracer())
    u0(x, y, z) = U0 * exp(z / max(d.Lz / 3, 1.0)) * tanh((y - 0.5d.Ly) / width)
    v0(x, y, z) = 0.05U0 * sin(2π * x / d.Lx) * exp(-((y - 0.5d.Ly)^2) / front_width^2)
    b0(x, y, z) = N² * z + B0 * tanh((y - 0.5d.Ly) / front_width) * exp(z / max(d.Lz / 2, 1.0))
    set!(model, u = u0, v = v0, b = b0)
    return model
end

"""
    dipole_vortex_builder(settings::RunSettings, arch)

Builds initial conditions modeled by two adjacent counter-rotating gaussian vortices to form a jet-propelling dipole vortex.
"""
function dipole_vortex_builder(settings::RunSettings, arch)
    d = settings.domain
    radius = 0.08 * min(d.Lx, d.Ly)
    amplitude = 1.2e10
    x1, y1 = 0.38d.Lx, 0.5d.Ly
    x2, y2 = 0.62d.Lx, 0.5d.Ly
    model = periodic_model(settings, arch; ν = 12.0)
    function u0(x, y, z)
        u1, _, _ = gaussian_vortex_velocity(x, y, x1, y1, radius, amplitude)
        u2, _, _ = gaussian_vortex_velocity(x, y, x2, y2, radius, -amplitude)
        return u1 + u2
    end
    function v0(x, y, z)
        _, v1, _ = gaussian_vortex_velocity(x, y, x1, y1, radius, amplitude)
        _, v2, _ = gaussian_vortex_velocity(x, y, x2, y2, radius, -amplitude)
        return v1 + v2
    end
    set!(model, u = u0, v = v0)
    return model
end

"""
    double_gyre_builder(settings::RunSettings, arch)

Builds initial conditions representing wind-driven double gyre streamfunctions common in generalized circulation simulations.
"""
function double_gyre_builder(settings::RunSettings, arch)
    d = settings.domain
    amplitude = 1.0e5
    model = periodic_model(settings, arch; ν = 16.0)
    u0(x, y, z) = -2π * amplitude / d.Ly * cos(2π * y / d.Ly) * sin(π * x / d.Lx)
    v0(x, y, z) =  π * amplitude / d.Lx * cos(π * x / d.Lx) * sin(2π * y / d.Ly)
    set!(model, u = u0, v = v0)
    return model
end

"""
    zonal_jet_instability_builder(settings::RunSettings, arch)

Configures an unstable Bickley jet perturbed with sinusoidal meanders to exhibit instability vortex formations.
"""
function zonal_jet_instability_builder(settings::RunSettings, arch)
    d = settings.domain
    width = 0.06 * d.Ly
    U0 = 1.1
    meander = 0.08 * d.Ly
    model = periodic_model(settings, arch; ν = 8.0)
    yj(x) = 0.5d.Ly + meander * sin(4π * x / d.Lx)
    u0(x, y, z) = U0 * sech((y - yj(x)) / width)^2
    v0(x, y, z) = 0.12U0 * sin(6π * x / d.Lx) * exp(-((y - 0.5d.Ly)^2) / (2width^2))
    set!(model, u = u0, v = v0)
    return model
end

"""
    vortex_merger_builder(settings::RunSettings, arch)

Constructs initial conditions with two identical co-rotating vortices that naturally fold and merge.
"""
function vortex_merger_builder(settings::RunSettings, arch)
    d = settings.domain
    radius = 0.09 * min(d.Lx, d.Ly)
    amplitude = 1.0e10
    x1, y1 = 0.43d.Lx, 0.5d.Ly
    x2, y2 = 0.57d.Lx, 0.5d.Ly
    model = periodic_model(settings, arch; ν = 10.0)
    function u0(x, y, z)
        u1, _, _ = gaussian_vortex_velocity(x, y, x1, y1, radius, amplitude)
        u2, _, _ = gaussian_vortex_velocity(x, y, x2, y2, radius, amplitude)
        return u1 + u2
    end
    function v0(x, y, z)
        _, v1, _ = gaussian_vortex_velocity(x, y, x1, y1, radius, amplitude)
        _, v2, _ = gaussian_vortex_velocity(x, y, x2, y2, radius, amplitude)
        return v1 + v2
    end
    set!(model, u = u0, v = v0)
    return model
end

"""
    isolated_eddy_builder(settings::RunSettings, arch)

Generates an isolated, drifting stable eddy traveling within a minor uniform background drift.
"""
function isolated_eddy_builder(settings::RunSettings, arch)
    d = settings.domain
    radius = 0.09 * min(d.Lx, d.Ly)
    amplitude = 1.0e10
    x0, y0 = 0.25d.Lx, 0.52d.Ly
    Ubg = 0.16
    model = periodic_model(settings, arch; ν = 8.0)
    function u0(x, y, z)
        u, _, _ = gaussian_vortex_velocity(x, y, x0, y0, radius, amplitude)
        return Ubg + u
    end
    function v0(x, y, z)
        _, v, _ = gaussian_vortex_velocity(x, y, x0, y0, radius, amplitude)
        return v
    end
    set!(model, u = u0, v = v0)
    return model
end

const EXAMPLE_SPECS = Dict(
    "barotropic_jet" => ExampleSpec(
        "barotropic_jet",
        "Barotropic jet",
        "A depth-uniform eastward jet with a weak transverse perturbation that rolls up into meanders.",
        120.0,
        5 * 86_400.0,
        10,
        DomainSettings(Lx = 1_400e3, Ly = 900e3, Lz = 400.0, Nx = 128, Ny = 96, Nz = 12),
        :relative_vorticity,
        barotropic_jet_builder,
    ),
    "baroclinic_jet" => ExampleSpec(
        "baroclinic_jet",
        "Baroclinic jet",
        "A surface-intensified jet coupled to a buoyancy front so the dashboard can show baroclinic structure in real time.",
        90.0,
        4 * 86_400.0,
        8,
        DomainSettings(Lx = 1_200e3, Ly = 850e3, Lz = 600.0, Nx = 128, Ny = 96, Nz = 16),
        :buoyancy,
        baroclinic_jet_builder,
    ),
    "dipole_vortex" => ExampleSpec(
        "dipole_vortex",
        "Dipole vortex",
        "Two counter-rotating eddies that self-propel across the domain as a balanced vortex pair.",
        75.0,
        3 * 86_400.0,
        8,
        DomainSettings(Lx = 1_000e3, Ly = 800e3, Lz = 350.0, Nx = 120, Ny = 96, Nz = 10),
        :relative_vorticity,
        dipole_vortex_builder,
    ),
    "double_gyre" => ExampleSpec(
        "double_gyre",
        "Double gyre",
        "A canonical two-gyre circulation seeded from an analytic streamfunction.",
        120.0,
        5 * 86_400.0,
        10,
        DomainSettings(Lx = 1_400e3, Ly = 900e3, Lz = 400.0, Nx = 128, Ny = 96, Nz = 12),
        :relative_vorticity,
        double_gyre_builder,
    ),
    "zonal_jet_instability" => ExampleSpec(
        "zonal_jet_instability",
        "Zonal jet instability",
        "A narrow Bickley jet with an imposed meander that breaks down into vortices.",
        60.0,
        2.5 * 86_400.0,
        6,
        DomainSettings(Lx = 1_200e3, Ly = 700e3, Lz = 300.0, Nx = 144, Ny = 88, Nz = 10),
        :relative_vorticity,
        zonal_jet_instability_builder,
    ),
    "vortex_merger" => ExampleSpec(
        "vortex_merger",
        "Vortex merger",
        "Two co-rotating eddies initialized close enough to wrap around and merge.",
        75.0,
        2.5 * 86_400.0,
        6,
        DomainSettings(Lx = 900e3, Ly = 700e3, Lz = 300.0, Nx = 120, Ny = 96, Nz = 10),
        :relative_vorticity,
        vortex_merger_builder,
        ),
    "isolated_eddy" => ExampleSpec(
        "isolated_eddy",
        "Isolated eddy propagation",
        "A compact eddy embedded in a weak background current so its translation is obvious in the live view.",
        75.0,
        3 * 86_400.0,
        8,
        DomainSettings(Lx = 1_300e3, Ly = 800e3, Lz = 300.0, Nx = 128, Ny = 96, Nz = 10),
        :relative_vorticity,
        isolated_eddy_builder,
    ),
)

"""
    example_options()

Returns a dashboard-friendly array mapping string labels to keys for each available example specification configuration.
"""
function example_options()
    return [(label = spec.label, value = spec.id) for spec in values(EXAMPLE_SPECS)]
end

"""
    field_options(spec::ExampleSpec)

Selects and parses variable output options specifically validating whether the physics setup permits buoyancy plotting.
"""
function field_options(spec::ExampleSpec)
    return [option for option in DEFAULT_FIELD_OPTIONS if option.value != "buoyancy" || supports_buoyancy(spec)]
end

"""
    field_options()

Returns the complete standard list of available viewer options supported by the framework.
"""
field_options() = collect(DEFAULT_FIELD_OPTIONS)

"""
    field_options(example_id::AbstractString)

Retrieves field visualization option metadata specific to the string configuration ID.
"""
function field_options(example_id::AbstractString)
    return field_options(EXAMPLE_SPECS[String(example_id)])
end

"""
    backend_options()

Offers CPU and (when available) GPU backend configuration values for simulation compute selection.
"""
backend_options() = collect(DEFAULT_BACKEND_OPTIONS)
