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

## ----------------------------------------------------------------------------------------
## ---------------- Neural Ordinary Differential Equations -------------------
# Simulate dataset (with the solution as an array)
println("Configuring problem...")
Δt = 0.1
f(u, p, t) = u.^3*p 
u₀ = [2.0 0]
tspan = (0.0, 25.0)
t = tspan[1]:Δt:tspan[end]
typeof(tspan)
p = [-.1 -2.; 2. -.1] #This will be the TRUE solution
prob = ODEProblem(f, u₀, tspan, p)
println("Done!")
## Here we can define how to solve it https://diffeq.sciml.ai/stable/solvers/ode_solve/
println("Solving numerically....")
sol= solve(prob, Tsit5(), saveat = Δt) # Generate 'true data'
x = [x[1] for x in sol.u]
y = [x[2] for x in sol.u]
plot(sol.t, x, title="Simulated data", label='x')
plot!(sol.t, y, title="Simulated data", c=:red, label='y')
println("Done!")
## ============= Simplest ODE dx/dt = ax --> sol: e^a =========================
# function ode_example(u, α, t)
#     α*u
# end
ode_example(u, p, t) = p*t
println("Setting Parameters")
u₀ = 0.5
tspan = (0.0, 1.0)
p = 1.1 #This will be the TRUE solution
println("Simulating data (true solution)...")
datasize = 30 
tsteps = range(tspan[1], tspan[2], length = datasize)
prob_trueode = ODEProblem(ode_example, u₀, tspan, p)

## ============= DELETE USED FOR PRESENTATION =========================
ode_data = solve(prob_trueode, Tsit5(), saveat = tsteps) # Generate 'true data'
scatter(tsteps, ode_data', title="Simulated data")
    

## Build our neural network W₂*thanh(W₁x + b₁) 
# dudt2 = FastChain((x, p) -> x.^3, # First thing you will receive x and and will cube the input
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

# -------------------------
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

## ============

## #################### (its seems it is just for speed) Sensitivity algorithms (2:55:23) #########################
using DiffEqSensitivity

# Here we follow the steps for inference of parameters
function loss(p; salg = InterpolatingAdjoint()) # This is our loss function
    tmp_prob = remake(prob, p=p)
    tmp_sol = solve(tmp_prob, Tsit5(), saveat = tsteps, sensealg = salg)
    sum(abs2, Array(tmp_sol) - ode_data)
end

@time res = DiffEqFlux.sciml_train(p->loss(p; salg=BacksolveAdjoint()), pinit, ADAM(0.05), 
                    maxiters = 100,
                    cb = neuralode_callback)

