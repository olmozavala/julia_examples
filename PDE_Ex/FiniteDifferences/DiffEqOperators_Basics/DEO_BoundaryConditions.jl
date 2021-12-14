# https://github.com/SciML/DiffEqOperators.jl 
using DiffEqOperators
using Plots
using LinearAlgebra
using BandedMatrices, BlockBandedMatrices, SparseArrays
gr()

## Define our test space (internal nodes only)
Δx = 0.1
lastx = 3*π/2
x = Δx:Δx:lastx-Δx
u = sin.(x)
l = @layout[a b]

Δx2d = 0.7
x2d = Δx2d:Δx2d:lastx-Δx2d
f(x,y) = sin.(x) + sin.(y)
u2d = zeros(length(x2d), length(x2d))
for i = 1:length(x2d), j = 1:length(x2d)
    u2d[i,j] = f(x2d[i], x2d[j])
end
a = scatter(x, u, ms=5, title="Example values of u  \n boundary conditions  \n will be tested")
b = heatmap(u2d)
plot(a,b, layout=l, )

## =========== Boundary Conditions ========================
## =========== 1D                  ========================

## -------- Dirichlet with 0 at the boundaries
bc = Dirichlet0BC(Float64) # Defines Dirichlet BC with 0

scatter(x, u, ms=5)
scatter!([0 x... lastx]', bc*u, ms=3, c=:red, title="Dirichlet 0", legend=nothing)
## --------- Dirichlet -------
bc = DirichletBC(-.5, .5) # Defines Dirichlet BC 

scatter(x, u, ms=5)
scatter!([0 x... lastx]', bc*u, ms=3, c=:red, title="Dirichlet l=.1, r=.8", legend=nothing)

## --------- Periodic BC -------
bc = PeriodicBC(Float64) # Defines Dirichlet BC with 0

scatter(x, u, ms=5)
scatter!([0 x... lastx]', bc*u, ms=3, c=:red, title="Periodic BC", legend=nothing)

## --------- Neumann  0-------
# bc = Neumann0BC(0.1) # Defines Neumann ∂f(a)=0 with 0 ath the boundaries. The only parameter is Δx
bc = Neumann0BC(0.1, 10) # Parameter is Δx and number of points

scatter(x, u, ms=5)
scatter!([0 x... lastx]', bc*u, ms=3, c=:red, title="Neumann0", legend=nothing)

## --------- Neumann 1-------
# Remember that the in the left side the value of derivative is substracted

# Defines Neumann with values of the derivatives in each side (on the left is always the negative of the value entered)
ldx = 2.0
rdx = 1.0
bc = NeumannBC((ldx, rdx), 0.1) 

scatter(x, u, ms=5)
scatter!([0 x... lastx]', bc*u, ms=3, c=:red, title="Neumann $(-ldx) in x[0] and $rdx in x[end]", legend=nothing)

## --------- Robin (function and derivative on the boundary)-------
# For each 'side' solve for u and u' this equation  use au(x) + b ∂u/∂n = c  
#The function receives a,b,c for each boundary
# bc = RobinBC((1.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.1) # Defines Dirichlet BC with 0` → 1u + 0 = 0 →  u(x) = 0
bc = RobinBC((1.0, 0.0, 2.0), (1.0, 0.0, 2.0), 0.1) # Defines Dirichlet BC with 2` → 1u + 0 = 2 →  u(x) = 2
# bc = RobinBC((0.0, 1.0, 1.0), (0.0, 1.0, 1.0), 0.1) # Defines NeumannBC with 1`

scatter(x, u, ms=5)
scatter!([0 x... lastx]', bc*u, ms=3, c=:red, title="Robin BC", legend=nothing)


## =========== 2D                  ========================
## --------- Dirichlet 2D  T, (Δx, Δy), order, (count_x, count_y) -------
dirx = 1.0
diry = 1.0
bcx = RobinBC{1}((1.0, 0.0, dirx), (1.0, 0.0, diry), (Δx2d), 1, size(u2d)) # Defines Dirichlet BC with 0`
dirx = 2.0
diry = 2.0
bcy = RobinBC{2}((1.0, 0.0, dirx), (1.0, 0.0, diry), (Δx2d), 1, size(u2d)) # Defines Dirichlet BC with 0`
bccomp = compose(bcx,bcy)

l = @layout[a b; c d]
a = heatmap(u2d, title="Original")
b = heatmap(bccomp*u2d, title="Both sides")
c = heatmap(bcx*u2d, title="Dirichlet 1 in x")
d = heatmap(bcy*u2d, title="Dirichlet 2 in y")
plot(a, b, c, d, layout=l, size=(800,800))

## --------- Neumann0 2D  T, (Δx, Δy), order, (count_x, count_y) -------
bcx, bcy = Neumann0BC(Float64, (.1, .1), 1, size(u2d)) # Defines Neumann0BC at two dimensions
bc = compose(bcx, bcy)

l = @layout[a b]
a = heatmap(u2d, title="Original")
b = heatmap(bc*u2d, title="Newmann 0 both sides")
plot(a, b, layout=l)

