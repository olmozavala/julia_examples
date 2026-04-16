# Plots Examples
# https://docs.juliaplots.org/latest/

# %% Basic setup
using Plots
using Random
using Distributions

# %% Making some simple plots with Plots
x = -2π:.1:3π
y = sin.(x)

# Basic line plot
plot(x, y, 
    title = "Simple Sine Wave",
    xlabel = "x",
    ylabel = "sin(x)",
    linewidth = 2,
    color = :blue)

# Scatter plot
scatter(x, y, 
    title = "Scatter Plot",
    xlabel = "x",
    ylabel = "sin(x)",
    markersize = 3,
    color = :red)

# Bar plot
bar(x[1:10:end], y[1:10:end], 
    title = "Bar Plot",
    xlabel = "x",
    ylabel = "sin(x)",
    color = :green)

# %% Multiple plots together
plot(x, [sin.(x), cos.(x), sin.(x).*cos.(x)], 
    title = "Multiple Functions",
    xlabel = "x",
    ylabel = "y",
    label = ["sin(x)" "cos(x)" "sin(x)×cos(x)"],
    linewidth = 2)

# %% Subplots
p1 = plot(x, sin.(x), title = "sin(x)", color = :blue)
p2 = plot(x, cos.(x), title = "cos(x)", color = :red)
p3 = plot(x, tan.(x), title = "tan(x)", color = :green)
p4 = plot(x, exp.(-x.^2), title = "Gaussian", color = :purple)

plot(p1, p2, p3, p4, layout = (2, 2), size = (800, 600))

# %% Heatmap
A = rand(50, 50)
heatmap(A, 
    title = "Random Matrix",
    colormap = :grays)

# %% 3D plots
x_3d = -2:0.1:2
y_3d = -2:0.1:2
z_3d = [exp(-(x^2 + y^2)) for x in x_3d, y in y_3d]

surface(x_3d, y_3d, z_3d, 
    title = "3D Surface Plot",
    colormap = :viridis)

# %% Scatter 3D
n_points = 100
x_scatter = randn(n_points)
y_scatter = randn(n_points)
z_scatter = randn(n_points)

scatter3d(x_scatter, y_scatter, z_scatter,
    title = "3D Scatter Plot",
    markersize = 3,
    color = 1:n_points,
    colormap = :viridis)

# %% ---- Recipes -----
@recipe function f(x)
    y = sin(x)
    seriestype --> :path
    x, y
end

plot(-2π:2π, title = "Custom Recipe")

# %% ---- Recipes with Distributions -----
function default_range(dist::Distribution, n = 4)
    μ, σ = mean(dist), std(dist)
    range(μ - n*σ, stop=μ + n*σ, length=100)
end

dist = Normal(10, 50)

@recipe function f(dist::Distribution, x = default_range(dist))
    y = map(xi -> pdf(dist,xi), x)
    seriestype --> :path
    x, y
end

plot(dist, title = "Normal Distribution PDF")

# %% ---- Animations -----
# Simple animation
anim = @animate for i in 1:100
    t = i * 0.1
    plot(x, sin.(x .+ t), 
        title = "Animated Sine Wave",
        ylim = (-1.5, 1.5),
        color = :blue,
        linewidth = 2)
end

gif(anim, "sine_animation.gif", fps = 10)

# %% ================= Lorenz Attractor ==============
println("Making Lorenz plot....")

Base.@kwdef mutable struct Lorenz
    dt::Float64 = 0.02
    σ::Float64 = 10
    ρ::Float64 = 28
    β::Float64 = 8/3
    x::Float64 = 1
    y::Float64 = 1
    z::Float64 = 1
end

function step!(l::Lorenz)
    dx = l.σ * (l.y - l.x);         l.x += l.dt * dx
    dy = l.x * (l.ρ - l.z) - l.y;   l.y += l.dt * dy
    dz = l.x * l.y - l.β * l.z;     l.z += l.dt * dz
end

attractor = Lorenz()

# Initialize a 3D plot with 1 empty series
plt = plot3d(
    1,
    xlim = (-30, 30),
    ylim = (-30, 30),
    zlim = (0, 60),
    title = "Lorenz Attractor",
    marker = 2,
)

# Build an animated gif by pushing new points to the plot, saving every 10th frame
@gif for i=1:1500
    step!(attractor)
    push!(plt, attractor.x, attractor.y, attractor.z)
end every 10

println("Done!") 