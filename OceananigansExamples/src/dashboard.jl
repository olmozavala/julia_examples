const FONT_URLS = [
    "https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=Space+Grotesk:wght@400;500;700&display=swap",
]

const CARD_STYLE = Dict(
    "background" => "rgba(9, 20, 31, 0.78)",
    "border" => "1px solid rgba(116, 186, 255, 0.14)",
    "borderRadius" => "20px",
    "boxShadow" => "0 24px 60px rgba(0, 0, 0, 0.24)",
)

const DEFAULT_OPENED_VARIABLES = [
    "depth_mean_relative_vorticity",
    "surface_buoyancy",
    "depth_mean_speed",
    "b",
    "u",
    "v",
    "w",
]

"""
    format_hours(seconds)

Converts the given `seconds` into an hours formatted string.
"""
format_hours(seconds) = @sprintf("%.2f h", seconds / 3600)
"""
    format_days(seconds)

Converts the given `seconds` into a days formatted string.
"""
format_days(seconds) = @sprintf("%.2f d", seconds / 86_400)

"""
    parse_float(x, fallback)

Safely attempts to parse `x` as a Float64, returning `fallback` upon any failure.
"""
function parse_float(x, fallback)
    x === nothing && return fallback
    x isa Number && return Float64(x)
    try
        return parse(Float64, string(x))
    catch
        return fallback
    end
end

"""
    parse_int(x, fallback)

Safely attempts to parse `x` as an Int, returning `fallback` upon any failure.
"""
function parse_int(x, fallback)
    x === nothing && return fallback
    x isa Integer && return Int(x)
    x isa Number && return round(Int, x)
    try
        return parse(Int, string(x))
    catch
        return fallback
    end
end

"""
    empty_figure(title, subtitle)

Constructs an empty transparent Plotly figure layout serving as a placeholder or fallback.
"""
function empty_figure(title, subtitle)
    return Dict(
        "data" => Any[],
        "layout" => Dict(
            "title" => Dict("text" => "<b>$title</b><br><span style='font-size:12px;color:#9ab7d5'>$subtitle</span>"),
            "paper_bgcolor" => "rgba(0,0,0,0)",
            "plot_bgcolor" => "rgba(7, 18, 29, 0.55)",
            "font" => Dict("family" => "Space Grotesk, sans-serif", "color" => "#e8f4ff"),
            "xaxis" => Dict("visible" => false),
            "yaxis" => Dict("visible" => false),
            "margin" => Dict("l" => 30, "r" => 30, "t" => 70, "b" => 30),
        ),
    )
end

"""
    dataset_attribute(ds, key, fallback="")

Extracts a global attribute `key` from the NCDataset `ds`, rendering `fallback` if not present.
"""
dataset_attribute(ds, key, fallback = "") = haskey(ds.attrib, key) ? ds.attrib[key] : fallback

"""
    run_label(path::String)

Builds a human-readable identifier from the attributes or dirname of a NetCDF `path`.
"""
function run_label(path::String)
    ds = NCDataset(path)
    try
        example_label = string(dataset_attribute(ds, "example_label", basename(dirname(path))))
        backend_label = string(dataset_attribute(ds, "backend_label", ""))
        run_name = basename(dirname(path))
        return isempty(backend_label) ? "$example_label | $run_name" : "$example_label | $backend_label | $run_name"
    finally
        close(ds)
    end
end

"""
    running_run_paths()

Retrieves a Set containing string paths to all currently active SimulationSessions.
"""
function running_run_paths()
    return Set(session.netcdf_path for session in list_sessions() if session.running)
end

"""
    readable_legacy_archive(path::String)

Checks whether an old or non-status NetCDF simulation archive possesses loadable variables like `time`.
"""
function readable_legacy_archive(path::String)
    try
        ds = NCDataset(path)
        try
            return haskey(ds, "time") && length(ds["time"]) > 0
        finally
            close(ds)
        end
    catch
        return false
    end
end

"""
    completed_archive(path::String)

Validates whether `path` serves a completed NetCDF simulation run directory by inspecting its locking markers.
"""
function completed_archive(path::String)
    isfile(path) || return false
    run_dir = dirname(path)
    isfile(run_marker_path(run_dir, COMPLETED_MARKER_FILENAME)) && return true
    isfile(run_marker_path(run_dir, RUNNING_MARKER_FILENAME)) && return false
    return readable_legacy_archive(path)
end

"""
    completed_run_paths()

Scans the global `RUNS_DIR` to collect file paths for all correctly fully resolved netCDF results.
"""
function completed_run_paths()
    mkpath(RUNS_DIR)
    blocked = running_run_paths()
    paths = String[]

    for run_dir in sort(readdir(RUNS_DIR; join = true); rev = true)
        path = joinpath(run_dir, NETCDF_FILENAME)
        if !(path in blocked) && completed_archive(path)
            push!(paths, path)
        end
    end

    return paths
end

"""
    run_file_options()

Builds an array of Dash dropdown options mapping human-readable labels to simulation netcdf file paths.
"""
function run_file_options()
    return [(label = run_label(path), value = path) for path in completed_run_paths()]
end

"""
    selected_run_value(paths, current_value)

Validates the user's `current_value` run choice against the existing simulated `paths`, updating selection dynamically if necessary.
"""
function selected_run_value(paths, current_value)
    isempty(paths) && return nothing
    current = current_value === nothing ? nothing : string(current_value)
    current in paths && return current
    return first(paths)
end

"""
    pretty_variable_label(name::String)

Converts underscore delimited `name` to spaced and capitalized legible titles.
"""
function pretty_variable_label(name::String)
    words = split(replace(name, '_' => ' '))
    return join(uppercasefirst.(words), " ")
end

"""
    is_plottable_variable(ds, name::String)

Confirms if a specified `name` key maps inside the NCDataset `ds` towards an array spanning x, y scales.
"""
function is_plottable_variable(ds, name::String)
    dims = Tuple(dimnames(ds[name]))
    return "time" in dims && any(startswith("x_"), dims) && any(startswith("y_"), dims)
end

"""
    field_variable_options(path::String)

Lists interactive plot options gathered safely sequentially from the plottable keys inside an NCDataset wrapper given by `path`.
"""
function field_variable_options(path::String)
    ds = NCDataset(path)
    try
        names = sort(filter(name -> is_plottable_variable(ds, name), String.(collect(keys(ds)))))
        return [(label = pretty_variable_label(name), value = name) for name in names]
    finally
        close(ds)
    end
end

"""
    choose_default_variable(options)

Suggests a default active trace based on the available `options` and fixed priority fallback targets.
"""
function choose_default_variable(options)
    values = [option.value for option in options]
    for candidate in DEFAULT_OPENED_VARIABLES
        candidate in values && return candidate
    end
    return isempty(values) ? nothing : first(values)
end

"""
    frame_marks(times)

Constructs dictionary representations dictating Dash slider intervals labeled by elapsed `times` tracking hours.
"""
function frame_marks(times)
    isempty(times) && return Dict(0 => "0 h")
    frame_count = length(times)
    step = max(1, ceil(Int, frame_count / 6))
    indices = unique(vcat(1:step:frame_count, frame_count))
    return Dict(index - 1 => format_hours(times[index]) for index in indices)
end

"""
    level_marks(z)

Constructs slider intervals detailing evenly-distributed vertical coordinates `z` across active volumes.
"""
function level_marks(z)
    isempty(z) && return Dict(0 => "2D")
    level_count = length(z)
    step = max(1, ceil(Int, level_count / 5))
    indices = unique(vcat(1:step:level_count, level_count))
    return Dict(index - 1 => @sprintf("%.0f m", z[index]) for index in indices)
end

"""
    variable_dims(ds, variable::String)

Reads string identifiers for multi-directional spatial spans from the NCDataset.
"""
function variable_dims(ds, variable::String)
    dims = Tuple(dimnames(ds[variable]))
    xdim = findfirst(dim -> startswith(dim, "x_"), dims)
    ydim = findfirst(dim -> startswith(dim, "y_"), dims)
    zdim = findfirst(dim -> startswith(dim, "z_"), dims)
    return (
        dims = dims,
        xdim = xdim === nothing ? nothing : dims[xdim],
        ydim = ydim === nothing ? nothing : dims[ydim],
        zdim = zdim === nothing ? nothing : dims[zdim],
    )
end

"""
    run_variable_summary(path::String, variable::String)

Extrapolates and summarizes core parameters (z depth array, explicit dims, and active sequence timestamps) associated with a trace variable wrapper.
"""
function run_variable_summary(path::String, variable::String)
    ds = NCDataset(path)
    try
        info = variable_dims(ds, variable)
        times = Float64.(Array(ds["time"][:]))
        z = info.zdim === nothing ? Float64[] : Float64.(Array(ds[info.zdim][:]))
        return (times = times, z = z, info = info)
    finally
        close(ds)
    end
end

"""
    frame_index(frame_value, frame_count)

Secures safely clamped bounded timeline targets tracking simulation steps.
"""
function frame_index(frame_value, frame_count)
    frame_count == 0 && return 0
    return clamp(parse_int(frame_value, frame_count - 1), 0, frame_count - 1)
end

"""
    level_index(level_value, level_count)

Secures bounded clamping targeting vertical elevation cross-sections across simulated nodes.
"""
function level_index(level_value, level_count)
    level_count == 0 && return 0
    return clamp(parse_int(level_value, level_count - 1), 0, level_count - 1)
end

"""
    variable_colors(name::String)

Identifies standard predefined thematic plotly colormaps and signed divergence hints suited specific variables.
"""
function variable_colors(name::String)
    lowercase_name = lowercase(name)
    if occursin("speed", lowercase_name)
        return ("Viridis", false)
    elseif occursin("vorticity", lowercase_name) || lowercase_name in ("u", "v", "w", "b", "surface_buoyancy")
        return ("RdBu", true)
    end
    return ("Thermal", false)
end

"""
    read_field_slice(path::String, variable::String, frame_value, level_value)

Extracts an optimized rendering segment block corresponding to specified variable arrays over isolated timelines slices.
"""
function read_field_slice(path::String, variable::String, frame_value, level_value)
    ds = NCDataset(path)
    try
        info = variable_dims(ds, variable)
        times = Float64.(Array(ds["time"][:]))
        t_index = frame_index(frame_value, length(times)) + 1
        z = info.zdim === nothing ? Float64[] : Float64.(Array(ds[info.zdim][:]))
        k_index = level_index(level_value, length(z)) + 1

        indices = Any[]
        for dim in info.dims
            if dim == info.zdim
                push!(indices, k_index)
            elseif dim == "time"
                push!(indices, t_index)
            else
                push!(indices, Colon())
            end
        end

        field = Array(ds[variable][indices...])
        x = Float64.(Array(ds[info.xdim][:])) ./ 1e3
        y = Float64.(Array(ds[info.ydim][:])) ./ 1e3
        z_value = isempty(z) ? nothing : z[k_index]

        return (
            x = x,
            y = y,
            field = permutedims(field),
            time = times[t_index],
            z_value = z_value,
            info = info,
            frame_number = t_index,
            frame_count = length(times),
        )
    finally
        close(ds)
    end
end

"""
    viewer_figure(path::String, variable::String, frame_value, level_value)

Packages parsed field traces, color properties, labels, and geometry shapes producing animated mapping figures suitable for display representations.
"""
function viewer_figure(path::String, variable::String, frame_value, level_value)
    !isfile(path) && return empty_figure("Completed Run Viewer", "Open a completed NetCDF file to animate saved fields.")

    slice = read_field_slice(path, variable, frame_value, level_value)
    colorscale, signed = variable_colors(variable)

    heatmap = Dict(
        "type" => "heatmap",
        "x" => slice.x,
        "y" => slice.y,
        "z" => slice.field,
        "colorscale" => colorscale,
        "zsmooth" => "best",
        "colorbar" => Dict("title" => Dict("text" => pretty_variable_label(variable))),
        "hovertemplate" => "x = %{x:.1f} km<br>y = %{y:.1f} km<br>value = %{z:.3e}<extra></extra>",
    )
    signed && (heatmap["zmid"] = 0.0)

    title_suffix = isnothing(slice.z_value) ? "" : @sprintf(" | z = %.0f m", slice.z_value)
    title = "<b>$(run_label(path))</b><br><span style='font-size:12px;color:#9ab7d5'>$(pretty_variable_label(variable)) | frame $(slice.frame_number) / $(slice.frame_count) | t = $(format_hours(slice.time))$(title_suffix)</span>"

    return Dict(
        "data" => Any[
            heatmap,
            Dict(
                "type" => "contour",
                "x" => slice.x,
                "y" => slice.y,
                "z" => slice.field,
                "showscale" => false,
                "hoverinfo" => "skip",
                "line" => Dict("width" => 0.45, "color" => "rgba(255,255,255,0.18)"),
                "contours" => Dict("coloring" => "none", "showlabels" => false),
            ),
        ],
        "layout" => Dict(
            "title" => Dict("text" => title),
            "paper_bgcolor" => "rgba(0,0,0,0)",
            "plot_bgcolor" => "rgba(7, 18, 29, 0.55)",
            "font" => Dict("family" => "Space Grotesk, sans-serif", "color" => "#e8f4ff"),
            "margin" => Dict("l" => 70, "r" => 28, "t" => 78, "b" => 55),
            "transition" => Dict("duration" => 250, "easing" => "cubic-in-out"),
            "uirevision" => "$(path)-$(variable)",
            "xaxis" => Dict("title" => "x (km)", "gridcolor" => "rgba(160, 199, 238, 0.10)", "zeroline" => false),
            "yaxis" => Dict("title" => "y (km)", "gridcolor" => "rgba(160, 199, 238, 0.10)", "zeroline" => false, "scaleanchor" => "x"),
        ),
    )
end

"""
    diagnostics_figure(path::String)

Queries and formats saved global metrics (like trace kinetic energy outputs over time) directly populating statistical summary lines mappings.
"""
function diagnostics_figure(path::String)
    !isfile(path) && return empty_figure("Saved diagnostics", "Open a completed NetCDF file to view run diagnostics.")

    ds = NCDataset(path)
    try
        haskey(ds, "kinetic_energy") || return empty_figure("Saved diagnostics", "This run does not contain scalar diagnostics.")
        time_days = Float64.(Array(ds["time"][:])) ./ 86_400
        traces = Any[
            Dict(
                "type" => "scatter",
                "mode" => "lines+markers",
                "name" => "Kinetic energy",
                "x" => time_days,
                "y" => Float64.(Array(ds["kinetic_energy"][:])),
                "line" => Dict("color" => "#50d4c8", "width" => 3),
                "marker" => Dict("size" => 6, "color" => "#ecfeff"),
            ),
        ]

        if haskey(ds, "enstrophy")
            push!(traces, Dict(
                "type" => "scatter",
                "mode" => "lines",
                "name" => "Enstrophy",
                "x" => time_days,
                "y" => Float64.(Array(ds["enstrophy"][:])),
                "yaxis" => "y2",
                "line" => Dict("color" => "#ffa24c", "width" => 2, "dash" => "dash"),
            ))
        end

        if haskey(ds, "max_speed")
            push!(traces, Dict(
                "type" => "scatter",
                "mode" => "lines",
                "name" => "Max speed",
                "x" => time_days,
                "y" => Float64.(Array(ds["max_speed"][:])),
                "line" => Dict("color" => "#f6d56f", "width" => 2),
            ))
        end

        return Dict(
            "data" => traces,
            "layout" => Dict(
                "title" => Dict("text" => "<b>Saved diagnostics</b><br><span style='font-size:12px;color:#9ab7d5'>Diagnostics reloaded from the completed NetCDF output.</span>"),
                "paper_bgcolor" => "rgba(0,0,0,0)",
                "plot_bgcolor" => "rgba(7, 18, 29, 0.55)",
                "font" => Dict("family" => "Space Grotesk, sans-serif", "color" => "#e8f4ff"),
                "margin" => Dict("l" => 65, "r" => 70, "t" => 78, "b" => 50),
                "legend" => Dict("orientation" => "h", "y" => 1.12, "x" => 0.0),
                "xaxis" => Dict("title" => "Simulation time (days)", "gridcolor" => "rgba(160, 199, 238, 0.10)"),
                "yaxis" => Dict("title" => "Kinetic energy / max speed", "gridcolor" => "rgba(160, 199, 238, 0.10)"),
                "yaxis2" => Dict("title" => "Enstrophy", "overlaying" => "y", "side" => "right"),
            ),
        )
    finally
        close(ds)
    end
end

"""
    metric_cards(sessions, completed_paths, opened_path)

Assembles HTML status component representations summarizing session lifecycle metric states.
"""
function metric_cards(sessions, completed_paths, opened_path)
    active_count = count(session -> session.running, sessions)
    failed_count = count(session -> session.failed, sessions)
    opened_label = isempty(opened_path) ? "None" : basename(dirname(opened_path))

    return [
        html_div(className = "metric-card", children = [
            html_div("Active runs", className = "metric-label"),
            html_div(string(active_count), className = "metric-value"),
        ]),
        html_div(className = "metric-card", children = [
            html_div("Completed NetCDFs", className = "metric-label"),
            html_div(string(length(completed_paths)), className = "metric-value"),
        ]),
        html_div(className = "metric-card", children = [
            html_div("Failed runs", className = "metric-label"),
            html_div(string(failed_count), className = "metric-value"),
        ]),
        html_div(className = "metric-card", children = [
            html_div("Opened run", className = "metric-label"),
            html_div(opened_label, className = "metric-value metric-value--small"),
        ]),
    ]
end

"""
    session_card(session::SimulationSession)

Generates individual log summaries displaying configuration details against a single logged application queue.
"""
function session_card(session::SimulationSession)
    timestamp = Dates.format(session.started_at, dateformat"yyyy-mm-dd HH:MM:SS")
    detail = session.failed ? session.error_message : session.status
    return html_div(className = "run-entry", children = [
        html_div(session.run_name, className = "run-entry-title"),
        html_div("$(replace(session.example_id, '_' => ' ')) | $(session.backend_label) | $timestamp", className = "run-entry-meta"),
        html_div(detail, className = "run-entry-copy"),
    ])
end

"""
    run_status_panel(sessions)

Returns formatted composite logs and active processes history mapping status messages corresponding to all executing simulation runs locally tracked.
"""
function run_status_panel(sessions)
    active = [session for session in sessions if session.running]
    recent = [session for session in sessions if !session.running][1:min(end, 3)]

    active_children = isempty(active) ? [html_div("No active model runs.", className = "run-entry-copy")] : session_card.(active)
    recent_children = isempty(recent) ? [html_div("No recently finished runs in this session.", className = "run-entry-copy")] : session_card.(recent)

    return [
        html_div("Run monitor", className = "status-chip"),
        html_h3("Background simulations", className = "status-title"),
        html_p("The dashboard only tracks run state while simulations execute. Visualization is loaded later from completed NetCDF files.", className = "status-copy"),
        html_div("Running now", className = "section-title section-title--small"),
        html_div(className = "run-list", children = active_children),
        html_div("Recent session history", className = "section-title section-title--small"),
        html_div(className = "run-list", children = recent_children),
        html_p("Run outputs are written under $(RUNS_DIR).", className = "status-copy status-copy--secondary"),
    ]
end

"""
    make_settings(example_id, dt, stop_time, output_interval, Lx, Ly, Lz, Nx, Ny, Nz)

Parses domain configurations cleanly producing strongly typed `RunSettings` properties.
"""
function make_settings(
    example_id,
    dt,
    stop_time,
    output_interval,
    Lx,
    Ly,
    Lz,
    Nx,
    Ny,
    Nz,
)
    spec = EXAMPLE_SPECS[string(example_id)]
    domain = DomainSettings(
        Lx = parse_float(Lx, spec.default_domain.Lx),
        Ly = parse_float(Ly, spec.default_domain.Ly),
        Lz = parse_float(Lz, spec.default_domain.Lz),
        Nx = parse_int(Nx, spec.default_domain.Nx),
        Ny = parse_int(Ny, spec.default_domain.Ny),
        Nz = parse_int(Nz, spec.default_domain.Nz),
    )
    return RunSettings(
        example_id = spec.id,
        dt = parse_float(dt, spec.default_dt),
        stop_time = parse_float(stop_time, spec.default_stop_time),
        output_interval = max(1, parse_int(output_interval, spec.default_output_interval)),
        domain = domain,
        display_variable = spec.default_display_variable,
        preferred_backend = :gpu,
    )
end

"""
    control_group(label, component)

Wraps an elementary interface element within categorized label containers securing structured application input layout patterns.
"""
function control_group(label, component)
    return html_div(className = "control-group", children = [
        html_label(label, className = "control-label"),
        component,
    ])
end

"""
    numeric_input(id, value; min=nothing, step="any")

Wraps a uniform numeric property component. 
"""
function numeric_input(id, value; min = nothing, step = "any")
    kwargs = Dict{Symbol, Any}(
        :id => id,
        :type => "number",
        :value => value,
        :step => step,
        :className => "control-input",
        :debounce => true,
    )
    min === nothing || (kwargs[:min] = min)
    return dcc_input(; kwargs...)
end

"""
    build_layout()

Structures the top-level HTML application canvas grouping control interfaces, status components and canvas graphing grids locally mapped.
"""
function build_layout()
    default_spec = EXAMPLE_SPECS["barotropic_jet"]
    d = default_spec.default_domain
    options = run_file_options()

    return html_div(className = "app-shell", children = [
        dcc_store(id = "launch-store"),
        dcc_store(id = "open-run-store", data = Dict("path" => "", "nonce" => 0)),
        dcc_interval(id = "refresh-interval", interval = 1_500, n_intervals = 0),
        dcc_interval(id = "animation-interval", interval = 700, n_intervals = 0),
        html_div(className = "hero-band", children = [
            html_div(className = "hero-copy", children = [
                html_div("Oceananigans / Dash / Julia", className = "eyebrow"),
                html_h1("Ocean Flow Lab", className = "hero-title"),
                html_p(
                    "Launch idealized ocean simulations in the background, save every run to NetCDF, and reopen completed outputs for post-run animation.",
                    className = "hero-description",
                ),
            ]),
            html_div(className = "hero-callout", children = [
                html_div("Run archive", className = "callout-label"),
                html_div("The main viewer no longer follows the live model. It opens only completed NetCDF files from the runs directory.", className = "callout-value"),
            ]),
        ]),
        html_div(className = "content-grid", children = [
            html_div(className = "controls-panel", style = CARD_STYLE, children = [
                html_h2("Launch run", className = "section-title"),
                control_group(
                    "Flow configuration",
                    dcc_dropdown(
                        id = "example-dropdown",
                        options = sort(example_options(); by = option -> option.value),
                        value = default_spec.id,
                        clearable = false,
                        className = "control-dropdown",
                    ),
                ),
                html_div(id = "example-description", className = "example-description", children = default_spec.description),
                html_div(className = "control-grid control-grid--two", children = [
                    control_group("dt (s)", numeric_input("dt-input", default_spec.default_dt; min = 1, step = 5)),
                    control_group("Stop time (s)", numeric_input("stop-time-input", default_spec.default_stop_time; min = 60, step = 300)),
                    control_group("Output cadence (steps)", numeric_input("output-interval-input", default_spec.default_output_interval; min = 1, step = 1)),
                ]),
                html_h3("Domain", className = "section-title section-title--small"),
                html_div(className = "control-grid control-grid--three", children = [
                    control_group("Lx (m)", numeric_input("Lx-input", d.Lx; min = 10_000, step = 10_000)),
                    control_group("Ly (m)", numeric_input("Ly-input", d.Ly; min = 10_000, step = 10_000)),
                    control_group("Lz (m)", numeric_input("Lz-input", d.Lz; min = 10, step = 10)),
                    control_group("Nx", numeric_input("Nx-input", d.Nx; min = 16, step = 8)),
                    control_group("Ny", numeric_input("Ny-input", d.Ny; min = 16, step = 8)),
                    control_group("Nz", numeric_input("Nz-input", d.Nz; min = 4, step = 2)),
                ]),
                html_button("Start run", id = "start-run-button", n_clicks = 0, className = "launch-button"),
            ]),
            html_div(className = "visuals-panel", children = [
                html_div(id = "metric-strip", className = "metric-grid", children = metric_cards(list_sessions(), completed_run_paths(), "")),
                html_div(className = "visual-card", style = CARD_STYLE, children = [
                    html_h2("Open completed run", className = "section-title"),
                    html_div(className = "viewer-control-row", children = [
                        control_group(
                            "NetCDF run",
                            dcc_dropdown(
                                id = "run-file-dropdown",
                                options = options,
                                value = isempty(options) ? nothing : first(options).value,
                                clearable = false,
                                className = "control-dropdown",
                            ),
                        ),
                        html_button("Open", id = "open-run-button", n_clicks = 0, className = "launch-button launch-button--secondary"),
                    ]),
                    html_div(className = "viewer-control-grid", children = [
                        control_group("Variable", dcc_dropdown(id = "viewer-variable-dropdown", options = Any[], value = nothing, clearable = false, className = "control-dropdown")),
                        control_group(
                            "Playback",
                            dcc_radioitems(
                                id = "timeline-mode",
                                options = [(label = "Pause", value = "pause"), (label = "Animate", value = "play")],
                                value = "pause",
                                className = "timeline-mode",
                                labelStyle = Dict("display" => "inline-flex", "alignItems" => "center", "gap" => "6px", "marginRight" => "14px"),
                                inputStyle = Dict("marginRight" => "6px"),
                            ),
                        ),
                    ]),
                    html_div(id = "opened-run-caption", className = "timeline-label", children = "Open a completed run to load variables and timesteps."),
                    control_group("Vertical level", dcc_slider(id = "level-slider", min = 0, max = 0, step = 1, value = 0, marks = Dict(0 => "2D"))),
                    html_div(id = "level-label", className = "timeline-label", children = "2D output"),
                    control_group("Timestep", dcc_slider(id = "frame-slider", min = 0, max = 0, step = 1, value = 0, marks = Dict(0 => "0 h"))),
                    html_div(id = "frame-label", className = "timeline-label", children = "No completed run opened."),
                    dcc_graph(id = "flow-graph", figure = empty_figure("Completed Run Viewer", "Open a completed NetCDF file to animate saved fields."), className = "flow-graph", animate = false, config = Dict("displayModeBar" => false)),
                ]),
                html_div(className = "visual-card", style = CARD_STYLE, children = [
                    dcc_graph(id = "diagnostics-graph", figure = empty_figure("Saved diagnostics", "Open a completed NetCDF file to view run diagnostics."), config = Dict("displayModeBar" => false)),
                ]),
            ]),
            html_div(className = "status-panel", style = CARD_STYLE, children = [
                html_div(id = "run-status", children = run_status_panel(list_sessions())),
            ]),
        ]),
    ])
end

"""
    build_app()

Assembles application components linking internal callback behaviors spanning inputs to background triggers wrapping the Dash dashboard architecture.
"""
function build_app()
    app = dash(external_stylesheets = FONT_URLS)
    app.layout = build_layout()

    callback!(
        app,
        Output("launch-store", "data"),
        Input("start-run-button", "n_clicks"),
        State("example-dropdown", "value"),
        State("dt-input", "value"),
        State("stop-time-input", "value"),
        State("output-interval-input", "value"),
        State("Lx-input", "value"),
        State("Ly-input", "value"),
        State("Lz-input", "value"),
        State("Nx-input", "value"),
        State("Ny-input", "value"),
        State("Nz-input", "value"),
    ) do n_clicks, example_id, dt, stop_time, output_interval, Lx, Ly, Lz, Nx, Ny, Nz
        (n_clicks isa Nothing || n_clicks < 1) && return Dict("status" => "idle")
        settings = make_settings(example_id, dt, stop_time, output_interval, Lx, Ly, Lz, Nx, Ny, Nz)
        session = start_run!(settings)
        return Dict("session_id" => session.id, "path" => session.netcdf_path, "status" => session.status)
    end

    callback!(
        app,
        Output("dt-input", "value"),
        Output("stop-time-input", "value"),
        Output("output-interval-input", "value"),
        Output("Lx-input", "value"),
        Output("Ly-input", "value"),
        Output("Lz-input", "value"),
        Output("Nx-input", "value"),
        Output("Ny-input", "value"),
        Output("Nz-input", "value"),
        Output("example-description", "children"),
        Input("example-dropdown", "value"),
    ) do example_id
        spec = EXAMPLE_SPECS[string(example_id)]
        d = spec.default_domain
        return (
            spec.default_dt,
            spec.default_stop_time,
            spec.default_output_interval,
            d.Lx,
            d.Ly,
            d.Lz,
            d.Nx,
            d.Ny,
            d.Nz,
            spec.description,
        )
    end

    callback!(
        app,
        Output("metric-strip", "children"),
        Output("run-status", "children"),
        Output("run-file-dropdown", "options"),
        Input("refresh-interval", "n_intervals"),
        Input("launch-store", "data"),
        State("open-run-store", "data"),
    ) do _, _, open_data
        sessions = list_sessions()
        paths = completed_run_paths()
        opened_path = open_data isa AbstractDict ? string(haskey(open_data, "path") ? open_data["path"] : (haskey(open_data, :path) ? open_data[:path] : "")) : ""
        return metric_cards(sessions, paths, opened_path), run_status_panel(sessions), run_file_options()
    end

    callback!(
        app,
        Output("run-file-dropdown", "value"),
        Input("refresh-interval", "n_intervals"),
        Input("launch-store", "data"),
        State("run-file-dropdown", "options"),
        State("run-file-dropdown", "value"),
    ) do _, _, options, current_run_value
        values = String[]
        if options isa AbstractVector
            for option in options
                if option isa AbstractDict
                    val = haskey(option, "value") ? option["value"] : (haskey(option, :value) ? option[:value] : "")
                    push!(values, string(val))
                elseif option isa NamedTuple && hasproperty(option, :value)
                    push!(values, string(getproperty(option, :value)))
                end
            end
        end
        return selected_run_value(values, current_run_value)
    end

    callback!(
        app,
        Output("open-run-store", "data"),
        Input("open-run-button", "n_clicks"),
        State("run-file-dropdown", "value"),
    ) do n_clicks, selected_path
        (n_clicks isa Nothing || n_clicks < 1 || selected_path === nothing) && return Dict("path" => "", "nonce" => 0)
        return Dict("path" => string(selected_path), "nonce" => parse_int(n_clicks, 0))
    end

    callback!(
        app,
        Output("viewer-variable-dropdown", "options"),
        Output("viewer-variable-dropdown", "value"),
        Output("opened-run-caption", "children"),
        Input("open-run-store", "data"),
    ) do open_data
        path = open_data isa AbstractDict ? string(haskey(open_data, "path") ? open_data["path"] : (haskey(open_data, :path) ? open_data[:path] : "")) : ""
        if isempty(path) || !isfile(path)
            return Any[], nothing, "Open a completed run to load variables and timesteps."
        end

        options = field_variable_options(path)
        value = choose_default_variable(options)
        return options, value, "Opened $(run_label(path))"
    end

    callback!(
        app,
        Output("frame-slider", "max"),
        Output("frame-slider", "marks"),
        Output("level-slider", "max"),
        Output("level-slider", "marks"),
        Input("open-run-store", "data"),
        Input("viewer-variable-dropdown", "value"),
    ) do open_data, variable
        path = open_data isa AbstractDict ? string(haskey(open_data, "path") ? open_data["path"] : (haskey(open_data, :path) ? open_data[:path] : "")) : ""
        if isempty(path) || !isfile(path) || variable === nothing
            return 0, Dict(0 => "0 h"), 0, Dict(0 => "2D")
        end

        summary = run_variable_summary(path, string(variable))
        return (
            max(0, length(summary.times) - 1),
            frame_marks(summary.times),
            max(0, length(summary.z) - 1),
            level_marks(summary.z),
        )
    end

    callback!(
        app,
        Output("frame-slider", "value"),
        Output("level-slider", "value"),
        Input("open-run-store", "data"),
        Input("viewer-variable-dropdown", "value"),
        Input("animation-interval", "n_intervals"),
        Input("timeline-mode", "value"),
        State("frame-slider", "value"),
        State("frame-slider", "max"),
        State("level-slider", "value"),
    ) do _, _, _, mode, current_value, max_value, current_level
        triggered = isempty(callback_context().triggered) ? "" : callback_context().triggered[1].prop_id
        current = parse_int(current_value, 0)
        max_index = parse_int(max_value, 0)
        level = parse_int(current_level, 0)

        if triggered in ("open-run-store.data", "viewer-variable-dropdown.value")
            return 0, 0
        end

        mode == "play" || return current, level
        max_index == 0 && return 0, level
        return (current >= max_index ? 0 : current + 1), level
    end

    callback!(
        app,
        Output("level-label", "children"),
        Output("frame-label", "children"),
        Input("open-run-store", "data"),
        Input("viewer-variable-dropdown", "value"),
        Input("frame-slider", "value"),
        Input("level-slider", "value"),
    ) do open_data, variable, frame_value, level_value
        path = open_data isa AbstractDict ? string(haskey(open_data, "path") ? open_data["path"] : (haskey(open_data, :path) ? open_data[:path] : "")) : ""
        if isempty(path) || !isfile(path) || variable === nothing
            return "2D output", "No completed run opened."
        end

        summary = run_variable_summary(path, string(variable))
        frame_idx = frame_index(frame_value, length(summary.times)) + 1
        level_idx = level_index(level_value, length(summary.z)) + 1

        level_caption = isempty(summary.z) ? "2D output" : @sprintf("Vertical level %d / %d | z = %.0f m", level_idx, length(summary.z), summary.z[level_idx])
        frame_caption = isempty(summary.times) ? "No saved timesteps." : "Frame $(frame_idx) / $(length(summary.times)) | t = $(format_hours(summary.times[frame_idx]))"

        return level_caption, frame_caption
    end

    callback!(
        app,
        Output("flow-graph", "figure"),
        Output("diagnostics-graph", "figure"),
        Input("open-run-store", "data"),
        Input("viewer-variable-dropdown", "value"),
        Input("frame-slider", "value"),
        Input("level-slider", "value"),
    ) do open_data, variable, frame_value, level_value
        path = open_data isa AbstractDict ? string(haskey(open_data, "path") ? open_data["path"] : (haskey(open_data, :path) ? open_data[:path] : "")) : ""
        if isempty(path) || !isfile(path) || variable === nothing
            return (
                empty_figure("Completed Run Viewer", "Open a completed NetCDF file to animate saved fields."),
                empty_figure("Saved diagnostics", "Open a completed NetCDF file to view run diagnostics."),
            )
        end

        return viewer_figure(path, string(variable), frame_value, level_value), diagnostics_figure(path)
    end

    return app
end

"""
    run_app(; host="127.0.0.1", port=8050, debug=true, hot_reload=true, ...)

Fires off server sockets dispatching application interactions utilizing predefined TCP loops effectively serving dashboard endpoints.
"""
function run_app(;
    host = "127.0.0.1",
    port = 8050,
    debug = true,
    hot_reload = true,
    hot_reload_interval = 0.5,
    hot_reload_watch_interval = 0.5,
    hot_reload_max_retry = 60,
)
    app = build_app()
    run_server(
        app,
        host,
        port;
        debug,
        dev_tools_hot_reload = hot_reload,
        dev_tools_hot_reload_interval = hot_reload_interval,
        dev_tools_hot_reload_watch_interval = hot_reload_watch_interval,
        dev_tools_hot_reload_max_retry = hot_reload_max_retry,
    )
end
