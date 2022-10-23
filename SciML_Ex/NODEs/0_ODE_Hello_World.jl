# Usings
# https://diffeq.sciml.ai/stable/types/ode_types/
using DifferentialEquations
using Plots
gr()

## Making plot function
function plotsol(t, u₀, α, sol) 
     plot(t, t->u₀*exp(α*t), lw=5, label="True Solution!", xlabel="t", ylabel="u(t)")
     p1 = plot!(sol,lw=3,ls=:dash,title="Num sol",
          xaxis="Time (t)",yaxis="u(t) (in μm)",label="ODE approx!") # legend=false
     p2 = plot(sol.t, sol.u, xaxis = "u(t)", yaxis="T", label="u(t)")
     plt = plot(p1, p2, layout=@layout[a b])
     display(plt)
end

## ################ Simplest case  du/dt = f(t,u) with f(t,u) = αu#############
# **************** Solving normally *************
println("Configuring problem...")
α = 1.1
Δt = 0.1
f(u, p, t) = p*u 
u₀ = 1/2
tspan = (0.0, 1.0)
t = tspan[1]:Δt:tspan[end]
typeof(tspan)
prob = ODEProblem(f, u₀, tspan, α)
println("Done!")
## Here we can define how to solve it https://diffeq.sciml.ai/stable/solvers/ode_solve/
println("Solving numerically....")
sol = solve(prob,reltol=1e-6,saveat=Δt) # Specify the error tolerance (compare and see)
println("Done!")
##
println("Plotting...")
plotsol(t, u₀, α, sol)
println("Done!")

## ******* Solved with ModelingTlkit: Simplest case  du/dt = f(t,u) with f(t,u) = αu*************
using ModelingToolkit
using DifferentialEquations
using Plots

println("Defining ODE system...")
@parameters t, α
@variables u(t)

Dt = Differential(t)
eqs = [Dt(u) ~ α*u]

@named sys = ODESystem(eqs)
println("Done!")

## Iniital conditions and solving problem
println("Setting ODE Problem and solving...")
u₀ = 1/2
u_ic = [u => u₀]  # Here we map the states to the initial conditions
p = [α => 1.1]
tspan = (0.0, 1.0)
prob = ODEProblem(sys, u_ic, tspan, p)
sol = solve(prob,reltol=1e-6,saveat=Δt) # Specify the error tolerance (compare and see)
println("Done!")

##
println("ODE Solved! Plotting...")
t = sol.t
α = 1.1
# plotsol(t, u₀, 1.1, sol)
# plot(t, t -> u₀*exp(α*t), lw=5, label="True Solution!", xlabel="t", ylabel="u(t)")
plot(t, t -> u₀*exp(α*t), lw=5, label="True Solution!", xlabel="t", ylabel="u(t)")
p1 = plot!(sol,lw=3,ls=:dash,title="Num sol",
     xaxis="Time (t)",yaxis="u(t) (in μm)",label="ODE approx!") # legend=false
p2 = plot(sol.t, sol.u, xaxis = "u(t)", yaxis="T", label="u(t)")
plt = plot(p1, p2, layout=@layout[a b])
display(plt)

println("Done!")