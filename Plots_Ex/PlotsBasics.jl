## http://docs.juliaplots.org/latest/basics/
using Random
using Plots
using Colors
using Printf
gr()
# plotly()
# pyplot()
print("The current backend in use is: ", backend())

## Printing the possible attributes 
plotattr(:Plot)
plotattr(:Series)
plotattr(:Subplot)
plotattr(:Axis)
# Plot specific attribute information
plotattr("size")

# --------------- BASICS --------------
## Making some simple plots
x = -2π:.1:3π
y = sin.(x)
#(Ctrl+Shift+Enter for line by line)
plot(x,y)
scatter(x, y)
bar(x, y)
##
title = @sprintf "Title %0.2f" 0.341232
plot(x, [y cos.(x)], 
title = title, 
xaxis="x",
yaxis="y", 
label =["m1" "m2"],
lw = 3)

## --------------- Call by function--------------
plot(x, sin)

# --------------- SUBPLOTS --------------
## Subplots
l = @layout[a b]
p1 = plot(y)
p2 = plot(2*y)
plot(p1, p2, layout=l)
# or 
plot(p1, p2, layout=@layout[a b])
## Modifying sizes
plot(p1, p2, layout=@layout[a{.9w};b{.3h}])

## Save to file
plot(sin,x)
savefig("MyFig.png")

## ---- Animations -----
x = 0:.01:2*π
# println("Done!")
@gif for i in range(0, stop = 2π, length = 100)
    f(x) = sin.(x.+i)
    p = plot(x, f)
end


## ----- 3D Plot ---


## ----- Vector fields -------
quiver([1,2,3],[3,2,1],quiver=([1,1,1],[1,2,3])) 
dh = .1
x = -π/2:dh:π/2
X = (x * ones(length(x))')'
Y = x * ones(length(x))'
# Single gyre
U(x,y) = -sin.(y)
V(x,y) =  sin.(x) .* sin.(y./2 .+ π/2)
# 
l = @layout[a b c d]
p1 = heatmap(x,x, X, title="X", c = :delta)
p2 = heatmap(x,x,Y, title="Y", c = :delta)
p3 = heatmap(x,x,U(X,Y), title="U", c = :delta)
p4 = heatmap(x,x,V(X,Y), title="V", c = :delta)
plot(p1, p2, p3, p4, layout=l, colorbar=:none)
##

function plotvector(x,y,u,v,subsample)
    l = 1:length(vec(X))
    lt = shuffle(l)
    ls = lt[1:subsample:end]
    quiver(vec(X)[ls],vec(Y)[ls],quiver=(vec(U(x,y))[ls],vec(V(x,y))[ls]))
end

plotvector(X,Y,U,V,2)

## Attributes
plot(x, [sin.(x) .* 1.3, cos.(x)], 
        title="Two plots together $sin(α)", label=["l1" "l2"], 
        c = [:blue :green],
        lw = 2, xlabel="xlabel", ylabel="ylabel") 
# Add to previous plot
plot!(x, sin.(x).*1.1, c = :red, bc = :yellow, size = (400,300))

# ## ---- Images/Matrices Not Working ----
heatmap(randn(10,10))

## ---- Recipes -----
@recipe function f{Float64}(x)
    y = sin(x)
    seriestype --> :path  # there is always an attribute dictionary `d` available...
    x, y
end
plot(-2π:2π, f)

## ---- Recipes -----
using Distributions
pyplot(size=(400,250));
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