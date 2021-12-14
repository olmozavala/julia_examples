using DiffEqOperators
using Plots
using LinearAlgebra
using BandedMatrices, BlockBandedMatrices, SparseArrays

## Define our test space (internal nodes only)
Δx = 0.4
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
plot(a,b, layout=l)

## =========== Basic exmaple ========================
# Define your boundary condition 
bc = Dirichlet0BC(Float64) # Defines Dirichlet BC with 0
scatter(x, u, ms=5)
scatter!([0 x... lastx]', bc*u, ms=3, c=:red, title="Dirichlet 0", legend=nothing)

## Compute second order approximation to the second derivative
order = 2
deriv = 2

N = length(x)
A = CenteredDifference(deriv, order, Δx, N)

println("Size of A is: $(size(A)). With values:")
display(Matrix(A))
println("Applied matrix A with boundary conditions on current u:$(Array(A*bc*u))" )

## Second order approximation to the second derivative in 2D
Δxx = CenteredDifference{1}(2,2,Δx,N)
Δyy = CenteredDifference{2}(2,2,Δx,N)
Δ = Δxx + Δyy
print("Size of A2 is: $(size(Δ))")
display(Matrix(Δ))
display(Matrix(Δxx))
display(Matrix(Δyy))
heatmap(Matrix(Δ))

## Compute second order approximation to the second derivative UpwindDifference
order = 2
deriv = 1

N = length(x)
A = UpwindDifference(deriv, order, Δx, N, 1)

println("Size of A is: $(size(A)). With values:")
display(Matrix(A))
println("Applied matrix A with boundary conditions on current u:$(Array(A*bc*u))" )

## Second order approximation to the second derivative in 2D
Δxx = UpwindDifference{1}(2,2,Δx,N)
Δyy = UpwindDifference{2}(2,2,Δx,N)
Δ = Δxx + Δyy
print("Size of A2 is: $(size(Δ))")
display(Matrix(Δ))
display(Matrix(Δxx))
display(Matrix(Δyy))
heatmap(Matrix(Δ))