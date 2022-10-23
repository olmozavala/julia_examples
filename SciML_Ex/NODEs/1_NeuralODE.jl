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
using Flux

## ---------------- Neural Ordinary Differential Equations -------------------
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
datasize = 30
tsteps = range(tspan[1], tspan[2], length = datasize)

p = [1.2, 0.6, 0.3, 0.8] #This will be the TRUE solution
println("Simulating data (true solution)...")
prob_trueode = ODEProblem(lotka_volterra!, u₀, tspan, p)
ode_data = Array(solve(prob_trueode, Tsit5(), saveat = tsteps)) # Generate 'true data'
scatter(tsteps, ode_data', title="Simulated data")

## Build our neural network W₂*thanh(W₁x + b₁) 
dudt2 = FastChain((x, p) -> x,  
                    FastDense(2, 50, tanh),  # Dense layer
                    FastDense(50, 2)) # Dense layer without activation

# Build a a neural ode, whenever you want to calculate the derivative call the NN
neural_ode_f(u,p,t) = dudt2(u,p)  
pinit = initial_params(dudt2) # Initial parameters of your ODE come from the intial parameters of your NN
prob = ODEProblem(neural_ode_f, u₀, tspan, pinit)  # We solve a NN using numerical ODE solveh

## Here we follow the steps for inference of parameters
function loss(p) # This is our loss function
    tmp_prob = remake(prob, p=p)
    tmp_sol = solve(tmp_prob, Tsit5(), saveat = tsteps)
    sum(abs2, Array(tmp_sol) - ode_data)
end
println("Done")

# -----------Visualization only --------------
function neuralode_callback(p,l)
    @show l # Shows the lss
    tmp_prob = remake(prob, p=p)
    tmp_sol = solve(tmp_prob, Tsit5(), saveat = tsteps)
    fig = plot(tmp_sol)
    scatter!(fig, tsteps, ode_data')
    display(fig)
    false # A boolean tells you if the optimizaqtion continues or not (false means no, it won't stop)
end

loss(pinit) # just to print initial error (with 'random parameters')

# Starts with Adam then moves to BFGS
# Search for controlling the adjoints (https://diffeqflux.sciml.ai/dev/ControllingAdjoints/)
@time res = DiffEqFlux.sciml_train(loss, pinit, ADAM(0.05), maxiters = 100, cb = neuralode_callback)
@time res2 =  DiffEqFlux.sciml_train(loss, res.minimizer, BFGS(initial_stepnorm=0.01), cb = neuralode_callback, maxiters = 100)

println("Done!..")

## Plotting the last parameters found
res_prob = ODEProblem(neural_ode_f, u₀, tspan, res2.minimizer)  # We solve a NN using numerical ODE solveh
tmp_sol = solve(res_prob, Tsit5(), saveat = tsteps)
scatter(tsteps, ode_data', title="Simulated data")
plot!(tmp_sol)

## ------- Different initial conditions ---------
u₀ = [1.0, 1.2]
res_prob = ODEProblem(neural_ode_f, u₀, tspan, res2.minimizer)  # We solve a NN using numerical ODE solveh
tmp_sol = solve(res_prob, Tsit5(), saveat = tsteps)

prob_trueode = ODEProblem(lotka_volterra!, u₀, tspan, p)
ode_data = Array(solve(prob_trueode, Tsit5(), saveat = tsteps)) # Generate 'true data'
scatter(tsteps, ode_data', title="Simulated data")
plot!(tmp_sol)