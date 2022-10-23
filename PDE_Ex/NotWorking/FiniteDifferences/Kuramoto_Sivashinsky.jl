using ModelingToolkit, DiffEqOperators, DifferentialEquations
using Plots
using Random
using DiffEqOperators, OrdinaryDiffEq
# NOT WORKING YET
# https://benchmarks.sciml.ai/html/MOLPDE/ks_fdm_wpd.html

gr()
include("../../OZ_Tools/ExampleData.jl")

## ========= Kuramoto-Sivashinsky  uₜ + Δ²u + Δu + 1/2|∇u|²
n = 100# Number of spatial grid points
distance = 4*π
Δx = distance/n
x = range(Δx, step=Δx, length=n) 
u₀ = -cos.(x) .* (1.0 .+ sin.(x))
plot(u₀, title="u₀")

##
ord_deriv = 1
ord_approx = 5
Δ⁴ = CenteredDifference(4, ord_approx, Δx, n)
Δ² = CenteredDifference(2, ord_approx, Δx, n)
Δ = CenteredDifference(1, ord_approx, Δx, n)

# Boundary conditions  au + b ∂u/∂n = c  
# bc = RobinBC((1.0, 0.0, 0.0),(1.0, 0.0, 0.0), (Δx), 1, size(u₀)) # Direichlet = 0
# bc = Neumann0BC(Δx)
# bc = Dirichlet0BC(Float64)
bc = PeriodicBC(Float64)

scatter(x, u₀, label="Scalar field")
display(scatter!(transpose([x[1]-Δx x... x[end]+Δx]), bc*u₀, label="BC", markersize=2, markerstrokewidth=0, mcolor=:red, c=:red))

## Solving equation
step(u,p,t) = -(Δ⁴*bc*u - Δ²*bc*u -.5*Δ*bc*u.^2)
# step(u,p,t) = -(Δ⁴*bc*u + Δ²*bc*u + .5*(Δ*bc*u).*u)

t0 = 0.0
t1 = 1.0
Δt = .0001
println("Solving with Δt=$Δt and Δx=$Δx")
prob = ODEProblem(step, u₀, (t0, t1))
# sol = solve(prob,  RadauIIA5(autodiff=false), abstol=1e-14, reltol=1e-14);
sol = solve(prob,  RadauIIA5(autodiff=false), abstol=1e-4, reltol=1e-4);
println("Done!")

## Plot results
tsteps = Int(round(t1/Δt))
println("Plotting...")
for i=1:length(sol.t)
    display(plot(sol.u[i], title="Time $(sol.t[i])"))
    sleep(.1)
end
##
heatmap(x, sol.t, hcat(sol.u...)', xaxis="Space", yaxis="Time")