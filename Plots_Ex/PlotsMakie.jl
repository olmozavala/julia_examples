# https://github.com/JuliaPlots/Makie.jl
## http://docs.juliaplots.org/latest/basics/
using Random
using Makie
using Colors

## Making some simple plots
scene = Scene()
x = -2π:.1:3π
y = sin.(x)

scene = lines(x,y, color = :blue)
scene = scatter!(x,y, color = rand(length(x)))

## ---- Images/Matrices ----
x = sin.(-2*π:.1:2*π)
# y = x .* transpose(x)
A = rand(50, 50)
scene = heatmap(A, color = :greys)

## --- Subplots
using AbstractPlotting: hbox, vbox
using AbstractPlotting

data = rand(50, 100)
p1 = heatmap(data, interpolate = true)
p2 = heatmap(data, interpolate = false)
p3 = heatmap(data, interpolate = true)
p4 = heatmap(data, interpolate = false)

s = hbox( 
        vbox(
            title(p1, "interpolate = true";  textsize = 15),
            title(p2, "interpolate = false"; textsize = 15),
        ),
        vbox(
            title(p3, "interpolate = true";  textsize = 15),
            title(p4, "interpolate = false"; textsize = 15),
        ),
        vbox(
            title(p1, "interpolate = true";  textsize = 15),
            title(p2, "interpolate = false"; textsize = 15),
        ),
        vbox(
            title(p3, "interpolate = true";  textsize = 15),
            title(p4, "interpolate = false"; textsize = 15),
        ),
)
# p_scene = Scene() # Parent scene
# child_scene = Scene(p_scene)
# x = sin.(-2*π:.1:2*π)
# scene = bar(x,y, color = rand(length(x)))
# display(child_scene)

## Attributes
scene = lines(x, [sin.(x) .* 1.3, cos.(x)], 
        title="Two plots together $sin(α)", label=["l1" "l2"], 
        c = [:blue :green],
        lw = 2, xlabel="xlabel", ylabel="ylabel") 

display(scene)
# Add to previous plot
plot!(x, sin.(x).*1.1, c = :red, bc = :yellow, size = (400,300))


## ---- Recipes -----
using Plots
pyplot(size=(400,250));
@recipe function f(x)
    y = sin(x)
    seriestype --> :path  # there is always an attribute dictionary `d` available...
    x, y
end
plot(-2π:2π)

## ---- Recipes -----
using Plots
pyplot(size=(400,250));
using Distributions
function default_range(dist::Distribution, n = 4)
    μ, σ = mean(dist), std(dist)
    range(μ - n*σ, stop=μ + n*σ, length=100)
end
dist = Normal(10, 50)
# Build a recipe which acts on a custom type.
# Notice that the function apply_recipe is returned.
# The recipe macro is just a convenience to build apply_recipe definitions.
@recipe function f(dist::Distribution, x = default_range(dist))
    y = map(xi -> pdf(dist,xi), x)
    seriestype --> :path  # there is always an attribute dictionary `d` available...
                          # If the user didn't specify a seriestype, we choose :path
    x, y
end
# that was pretty easy!
plot(dist)

## ---- Animations -----


## ================= More complicated stuff ==============
# Examples: http://docs.juliaplots.org/latest/
# define the Lorenz attractor

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

# initialize a 3D plot with 1 empty series
plt = plot3d(
    1,
    xlim = (-30, 30),
    ylim = (-30, 30),
    zlim = (0, 60),
    title = "Lorenz Attractor",
    marker = 2,
)

# build an animated gif by pushing new points to the plot, saving every 10th frame
@gif for i=1:1500
    step!(attractor)
    push!(plt, attractor.x, attractor.y, attractor.z)
end every 10
println("Done!")
