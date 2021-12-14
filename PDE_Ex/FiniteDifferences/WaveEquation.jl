using DiffEqOperators, DifferentialEquations, OrdinaryDiffEq
using Plots
gr()
include("../../OZ_Tools/ExampleData.jl")

## ====== Solving 1D wave equation uₜₜ = c²uₓₓ 
# https://en.wikipedia.org/wiki/Wave_equation

Δx = .1
# x = Δx:Δx:2*π
x = 0:Δx:2*π
n = length(x)
u₀ = sin.(x)
# u₀ = getGaussShape(n, n/8)
du₀ = -sin.(x)
# du₀ = -ones(n)

# ------- Define Boundary conditions -------
bc = Dirichlet0BC(Float64)
# bc = PeriodicBC(Float64)
# Plots intial condition with boundary conditions
scatter(x, u₀, c=:blue, legend=nothing, title="u₀")
scatter!(x, du₀, c=:green, legend=nothing, title="du₀")
scatter!([x[1]-Δx x... x[end]+Δx]', bc*u₀, ms=2, c=:Red, title="u₀ with bc")

##
println("Defining the problem...")
order = 2 # approximation order
deriv = 2 # Order of derivative
c = .001

Δ = CenteredDifference{1}(deriv, order, Δx, n)

tspan = (0.0, 5.0)
Δt = .0001
tsteps = Int(tspan[2]/Δt)

println("Solving Wave equation uₜₜ = c²uₓₓ with Δt=$Δt and Δx=$Δx")
# Solve an ODE uₜₜ = c²uₓₓ

# --- Using second order ODE solver ----
function wave(ddu, du,u,p,t) 
    ddu = c^2*(Δ*bc*u) # z' = c²uₓₓ
end
prob = SecondOrderODEProblem(wave, du₀, u₀, tspan)

# --- Manually 'joining' u and du ----------------k
# function wave(du,u,p,t) 
#     # Because we need to add the variable z = du/dt, we need to add
#     # another verctor with the same size of the spatial dimension of U.
#     # In this case u' is on the first n elements, and u on the other n
#     du[1:n] = u[1:n] # z = u'
#     du[n+1:n*2] = c^2*(Δ*bc*u[n+1:n*2]) # z' = c²uₓₓ
# end
# v₀ = [du₀... u₀...]
# prob = ODEProblem(wave, v₀, tspan)

sol = solve(prob, saveat=Δt)
## Plot
for t = 1:Int(round(tsteps/10)):tsteps
    display(plot(x,sol.u[t][n+1:n*2], c=:blue, title="Solution time $(sol.t[t])", 
                legend=nothing, ylims=(-1.5,1.5)))
    sleep(.1)
end
# println("Done!")

## - integrator example
# integrator = init(prob, Tsit5(), saveat=Δt)
# for t in 1:10
#     display(plot(x, integrator.u[n+1:n*2], c=:blue, title="Solution time $(integrator.t))", legend=nothing))
#     step!(integrator)
# end
## ====== Solving 2D wave equation uₜₜ = c²uₓₓ 



# plot(sol[1,:,1], c=:blue, title="Solution", legend=nothing)
    #
# heatmap(sol[1,:,:]', title="Solution")
# yaxis!("T")
# xaxis!("x")
