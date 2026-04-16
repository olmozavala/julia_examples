# Makie Examples
# https://github.com/JuliaPlots/Makie.jl
# https://makie.juliaplots.org/stable/

# %% Basic setup
using Random
using CairoMakie
using Colors
using GeometryBasics

# %% -------------- Simple plots --------------
x = -2π:.1:3π
y = sin.(x)

# Create a figure and axis
fig = Figure()
ax = Axis(fig[1, 1], 
    title = "Simple Sine Wave",
    xlabel = "x",
    ylabel = "sin(x)")

# Add line plot
CairoMakie.lines!(ax, x, y, color = :blue, linewidth = 2)
# Add scatter plot with random colors 
colors = rand(length(x))  # Single values for colormap
CairoMakie.scatter!(ax, x, y, color = colors, markersize = 10, colormap = :viridis)
fig

# %% ----------- Scatter with different colors and sizes --------------
points = Point2f.(x, y .+ 1)
fig = Figure()
ax = Axis(fig[1, 1], title = "Scatter Plot Examples")
CairoMakie.scatter!(ax, points, color = 1:length(points), markersize = range(5, 15, length = length(points)),
    colormap = :thermal, label = "colored points")

CairoMakie.axislegend(ax)
fig

# %% ---------------- Scatterlines example (from latest docs) ----------------
fig = Figure()
ax = Axis(fig[1, 1], title = "Scatterlines Examples")

xs = LinRange(0, 10, 20)
ys = 0.5 .* sin.(xs)

CairoMakie.scatterlines!(ax, xs, ys, color = :red, label = "red line")
CairoMakie.scatterlines!(ax, xs, ys .- 1, color = xs, markercolor = :red, label = "colored line")
CairoMakie.scatterlines!(ax, xs, ys .- 2, markersize = LinRange(5, 15, 20), label = "varying size")
CairoMakie.scatterlines!(ax, xs, ys .- 3, marker = :cross, strokewidth = 1,
    markersize = 8, color = :orange, strokecolor = :black, label = "cross markers")

CairoMakie.axislegend(ax)
fig

# %% -------------- Heatmaps for images or matrices--------------
A = rand(50, 50)
fig = Figure()
ax = Axis(fig[1, 1], title = "Random Matrix Heatmap")
CairoMakie.heatmap!(ax, A, colormap = :viridis)
fig

# %% ---------------- Streamplots ----------------
struct FitzhughNagumo{T}
    ϵ::T
    s::T
    γ::T
    β::T
end

P = FitzhughNagumo(0.1, 0.0, 1.5, 0.8)

f(x, P::FitzhughNagumo) = Point2f(
    (x[1]-x[2]-x[1]^3+P.s)/P.ϵ,
    P.γ*x[1]-x[2] + P.β
)

f(x) = f(x, P)

fig = Figure()
ax = Axis(fig[1, 1], title = "Fitzhugh-Nagumo Streamplot")
streamplot!(ax, f, -1.5..1.5, -1.5..1.5, colormap = :magma)
fig

# %% ---------------- Subplots using modern Makie syntax ----------------
data = rand(50, 100)
fig = Figure()
ax1 = Axis(fig[1, 1], title = "interpolate = true")
ax2 = Axis(fig[1, 2], title = "interpolate = false")
ax3 = Axis(fig[2, 1], title = "interpolate = true")
ax4 = Axis(fig[2, 2], title = "interpolate = false")

CairoMakie.heatmap!(ax1, data, interpolate = true)
CairoMakie.heatmap!(ax2, data, interpolate = false)
CairoMakie.heatmap!(ax3, data, interpolate = true)
CairoMakie.heatmap!(ax4, data, interpolate = false)

fig

# %% ---------------- Multiple lines with legend ----------------
fig = Figure()
ax = Axis(fig[1, 1], 
    title = "Multiple Functions",
    xlabel = "x", 
    ylabel = "y")

CairoMakie.lines!(ax, x, sin.(x) .* 1.3, color = :blue, linewidth = 2, label = "1.3 × sin(x)")
CairoMakie.lines!(ax, x, cos.(x), color = :green, linewidth = 2, label = "cos(x)")
CairoMakie.lines!(ax, x, sin.(x).*1.1, color = :red, linewidth = 2, label = "1.1 × sin(x)")

CairoMakie.axislegend(ax)
fig

# %% ---------------- Bar plot example ----------------
fig = Figure()
ax = Axis(fig[1, 1], 
    title = "Bar Plot Example",
    xlabel = "Categories",
    ylabel = "Values")

categories = ["A", "B", "C", "D", "E"]
values = rand(5)
CairoMakie.barplot!(ax, 1:5, values, color = :steelblue)
ax.xticks = (1:5, categories)

fig

# %% ---------------- 3D Scatter plot ----------------
fig = Figure()
ax = Axis3(fig[1, 1], 
    title = "3D Scatter Plot",
    xlabel = "x",
    ylabel = "y", 
    zlabel = "z")

# Generate 3D data
n_points = 100
x_3d = randn(n_points)
y_3d = randn(n_points)
z_3d = randn(n_points)

# Create RGB colors explicitly to avoid RGBA conversion issues
colors = [RGB(rand(), rand(), rand()) for _ in 1:n_points]
CairoMakie.scatter!(ax, x_3d, y_3d, z_3d, color = colors, markersize = 10)
fig

# %% ---------------- Surface plot ----------------
fig = Figure()
ax = Axis3(fig[1, 1], 
    title = "Surface Plot",
    xlabel = "x",
    ylabel = "y", 
    zlabel = "z")

x_surf = -2:0.1:2
y_surf = -2:0.1:2
z_surf = [exp(-(x^2 + y^2)) for x in x_surf, y in y_surf]

CairoMakie.surface!(ax, x_surf, y_surf, z_surf, colormap = :viridis)
fig

# %% ---------------- Animation example (simple) ----------------
fig = Figure()
ax = Axis(fig[1, 1], 
    title = "Animated Plot",
    xlabel = "x",
    ylabel = "y",
    limits = (-2π, 2π, -1.5, 1.5))

# Create animation
record(fig, "animation.gif", 1:100) do frame
    t = frame * 0.1
    y_anim = sin.(x .+ t)
    
    # Clear previous plot
    empty!(ax)
    
    # Add new line
    CairoMakie.lines!(ax, x, y_anim, color = :blue, linewidth = 2)
end

println("Animation saved as 'animation.gif'")

# %% ---------------- Lorenz Attractor ----------------
println("Making Lorenz plot....")
Base.@kwdef mutable struct Lorenz
    dt::Float64 = 0.01 # Time step
    σ::Float64 = 10  # Prandtl number
    ρ::Float64 = 28 # Rayleigh number
    β::Float64 = 8/3 # Ratio of specific heats
    x::Float64 = 0.1 # Initial condition
    y::Float64 = 0.1 # Initial condition
    z::Float64 = 0.1 # Initial condition
end

function step!(l::Lorenz)
    dx = l.σ * (l.y - l.x);         l.x += l.dt * dx
    dy = l.x * (l.ρ - l.z) - l.y;   l.y += l.dt * dy
    dz = l.x * l.y - l.β * l.z;     l.z += l.dt * dz
end

# Generate Lorenz attractor data
attractor = Lorenz()
n_steps = 5000
trajectory = zeros(n_steps, 3)

for i in 1:n_steps
    trajectory[i, :] = [attractor.x, attractor.y, attractor.z]
    step!(attractor)
end

# Plot 3D trajectory
fig = Figure()
ax = Axis3(fig[1, 1], 
    title = "Lorenz Attractor",
    xlabel = "x",
    ylabel = "y", 
    zlabel = "z"
)

CairoMakie.lines!(ax, trajectory[:, 1], trajectory[:, 2], trajectory[:, 3], 
    color = 1:n_steps, colormap = :viridis, linewidth = 1)
println("Done!")
fig