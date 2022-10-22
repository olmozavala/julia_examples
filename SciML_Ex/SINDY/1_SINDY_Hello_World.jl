using DataDrivenDiffEq
using ModelingToolkit
using LinearAlgebra
using Plots
##

f(u) = u.^2 .+ 2.0u .- 1.0
x = randn(1, 100);
y = reduce(hcat, map(f, eachcol(x)));
scatter(x',y')

## Define basis
@variables u
basis = Basis(monomial_basis([u], 2), [u])
println(basis)

##
problem = DirectDataDrivenProblem(x, y)
sol = solve(problem, basis, STLSQ())
print("Done!")
##
println(sol)
println(result(sol))
println(sol.parameters)