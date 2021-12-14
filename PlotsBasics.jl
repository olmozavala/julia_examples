## http://docs.juliaplots.org/latest/basics/
using Random
using Plots
using Colors
gr()
# plotly()
print("The current backend in use is: ", backend())

## Printing the possible attributes 
plotattr(:Plot)
plotattr(:Series)
plotattr(:Subplot)
plotattr(:Axis)
# Plot specific attribute information
plotattr("size")

## Making some simple plots
x = -2π:.1:3π
y = sin.(x)

plot(x,y)
# Same as 
plot(sin,x)
scatter(x, y)
bar(x, y)

# Same as 

## Save to file
plot(sin,x)
savefig("MyFig.png")

## Attributes
plot(x, [sin.(x) .* 1.3, cos.(x)], 
        title="Two plots together $sin(α)", label=["l1" "l2"], 
        c = [:blue :green],
        lw = 2, xlabel="xlabel", ylabel="ylabel",
        ylims = (-2,2)) 
# Add to previous plot
plot!(x, sin.(x).*1.1, c = :red, bc = :yellow, size = (800,600))

## Subplots (run with single line, not sure why it doesnt work for run cell)
l = @layout [a b]
a = plot(sin.(x));
b = plot(cos.(x));
plot(a,b, layout=l)
a = heatmap(rand(10,10));
b = heatmap(rand(10,10));
plot(a,b, layout=@layout [a b])


# ## ---- Images/Matrices Not Working ----
heatmap(randn(10,10), clim=(0,1))

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