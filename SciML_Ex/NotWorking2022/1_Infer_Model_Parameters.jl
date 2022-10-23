# Usings
# https://live.juliacon.org/talk/C9FGPP   (Chris workshop)
# https://sciml.ai/
# https://github.com/SciML/SciMLTutorials.jl   (Tutorials)
# https://tutorials.sciml.ai/
# https://github.com/SciML/NeuralPDE.jl    (for PDE)
# https://diffeq.sciml.ai/dev/solvers/ode_solve/   (Here you can see the methods and suggests which one to use)
using DifferentialEquations
using Plots
using DiffEqFlux, Optim
using DiffEqSensitivity
using Flux
gr()

## ------------------ First example infer the model parameters (2:11:17 of Chris talk)-----------------------
# Simulate dataset (with the solution as an array)
function lotka_volterra!(du, u, p, t)
    r, w = u  # Obtain current values of rabbits and wolves from U
    α, β, γ, δ = p
    du[1] =  α*r - β*r*w # Prey
    du[2] =  γ*r*w - δ*w # Predator
end

println("Setting Parameters")
u₀ = [1.0, 1.0]
tspan = (0.0, 10.0)
p = [1.2, 0.6, 0.3, 0.8] #This will be the TRUE solution

println("Simulating data (true solution)...")
prob = ODEProblem(lotka_volterra!, u₀, tspan, p)
sol = solve(prob, saveat=0.1)  # Here you can choose your solver and the accuracy desired options are here: https://diffeq.sciml.ai/dev/basics/common_solver_opts/
# Here is where we get an array from the solution
dataset = Array(sol)

plot(sol)
scatter!(sol.t, Array(sol)')  # Transposed dataset, just because of the way plots and difeq saves the 
# Manual plotting
# scatter(sol.t, [x[1] for x in sol.u[:]], label="Prey")
# scatter!(sol.t, [x[2] for x in sol.u[:]], label="Predator", xaxis="Time", yaxis="Amount ")
##


println("Define a loss function")
## Define a loss function that it computes a new model with input paramteres and the temporal problem
function loss(p)
    tmp_prob = remake(prob, p=p)
    tmp_sol = solve(tmp_prob, Tsit5(), saveat=0.1)
    sum(abs2, Array(tmp_sol) - dataset), tmp_sol
end

# Define a plot callback so that we can see the intermediate results
# It is receiving the output forom loss
function plot_callback(p, l, tmp_sol)
    @show l
    tmp_prob = remake(prob, p=p)
    tmp_sol = solve(tmp_prob, saveat=0.1, verbose=false)
    fig = plot(tmp_sol, label="Approximated solution")
    scatter!(sol.t, dataset', label="True data")
    display(fig)
    false  # Don't return anything
end

println("Optimize parameters with BFGS...")
# Initial guess
pinit = [1.2, 0.8, 2.5, 0.8]
# res = DiffEqFlux.sciml_train(loss, pinit, ADAM(0.02), cb=plot_callback, maxiters=300)
res = DiffEqFlux.sciml_train(loss, pinit, ADAM(0.02), maxiters=300)
println("Parameters that it found: ", res.minimizer)

## Plot results 
prob = ODEProblem(lotka_volterra!, u₀, tspan, res.minimizer)
sol = solve(prob, saveat=0.1)  # Here you can choose your solver and the accuracy desired options are here: https://diffeq.sciml.ai/dev/basics/common_solver_opts/
scatter(sol.t, dataset', label="Data")
plot!(sol, label="Approx")