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

## ------------------ First example infer the model parameters (2:11:17 of Chris talk)-----------------------
# Simulate dataset (with the solution as an array)
function lotka_volterra!(du, u, p, t)
    r, w = u  # Obtain current values of rabbits and wolves from U
    α, β, γ, δ = p
    du[1] = dr = α*r - β*r*w
    du[2] = dw = γ*r*w - δ*w
end

print("Parameters")
u₀ = [1.0, 1.0]
tspan = (0.0, 10.0)
p = [1.2, 0.6, 0.3, 0.8] #This will be the TRUE solution

print("Simulating data (true solution)")
prob = ODEProblem(lotka_volterra!, u₀, tspan, p)
sol = solve(prob, saveat=0.1)  # Here you can choose your solver and the accuracy desired options are here: https://diffeq.sciml.ai/dev/basics/common_solver_opts/
# Here is where we get an array from the solution
dataset = Array(sol)
plot(sol)
scatter!(sol.t, dataset')  # Transposed dataset, just because of the way plots and difeq saves the 

print("Define a loss function")
## Define a loss function that it computes a new model with input paramteres and the temporal problem
function loss(p)
    tmp_prob = remake(prob, p=p)
    tmp_sol = solve(tmp_prob, saveat=0.1)
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

print("Optimize parameters with BFGS")
# Initial guess
pinit = [1.2, 0.8, 2.5, 0.8]
res = DiffEqFlux.sciml_train(loss, pinit, ADAM(0.02), cb=plot_callback, maxiters=300)
print("Parameters that it found: ", res.minimizer)

# Here is the same but using a global optimizer, will need to read the docs. 
# using BlackBoxOptim  
# res = DiffEqFlux.sciml_train(loss, pinit, DiffEqFlux.BBO(), lower_bounds=0.2ones(4), upper_bounds=4.0ones(4))


## ----------------------------------------------------------------------------------------
## --------------Bayessian inference (predict parameters with uncertainty) (2:23 from workshop) ------------------------------------------------
# Website: https://turing.ml/dev/tutorials/10-bayesiandiffeq/     
using DiffEqFlux, OrdinaryDiffEq, Flux, Optim, Plots, Turing, MCMCChains, StatsPlots
Turing.setadbackend(:forwarddiff)  # Define which type of differentitation we want

# Making Maximul Likelihood
# For OZ, example of @model: https://julialang.org/blog/2017/08/dsl/
@model function fitlv(data)
    # Define prior distributions for each parameter
    σ ~ InverseGamma(2, 3)
    α ~ truncated(Normal(1.3, 0.5), 0.5, 2.5)  
    β ~ truncated(Normal(1.2, 0.5), 0, 2)  
    γ ~ truncated(Normal(2.7, 0.5), 1, 4)
    δ ~ truncated(Normal(1.3, 0.5), 0, 2)

    # Then you define how to compute the model
    p = [α, β, γ, δ]
    prob = ODEProblem(lotka_volterra!, u₀, (0.0,10.0), p)  # Notice that we cant use tspan here because it is not defined inside the function
    predicted = solve(prob, Tsit5(), saveat=0.1)

    # Check predictions with the data
    for i = 1:length(predicted)
        # MvNormal receives Covariance matrix and mean
        data[:,i] ~ MvNormal(predicted[i], σ) #  
    end
end

# p = [1.2, 0.6, 0.3, 0.8] These are the true values
# (results are not very good, we need to reviewed it) not sure what is σ and it appears in the plots as parameter
model = fitlv(dataset) # We call our function to create the model
chain = sample(model, NUTS(.35), 3000) # 
plot(chain)


#########################
## ----------------------------------------------------------------------------------------
## ---------------- Neural Differential Equations 2:30 from talk  -------------------
# We are working with A*u^3 in 2 dimensions.
# Here we want to find the model for us

u0 = Float32[2.0; 0.0]
datasize = 30
tspan = (0.0f0, 1.5f0)
tsteps = range(tspan[1], tspan[2], length = datasize)

# Simulate data
function trueODEfunc(du, u, p, t)
    true_A = [-0.1 2.0; -2.0 -0.1]
    du .= ((u.^3)'true_A)'
end

prob_trueode = ODEProblem(trueODEfunc, u0, tspan)
ode_data = Array(solve(prob_trueode, Tsit5(), saveat = tsteps)) # Generate 'true data'
    
# Build our neural network W₂*thanh(W₁x + b₁) 
dudt2 = FastChain((x, p) -> x.^3, # First thing you will receive x and and will cube the input
                    FastDense(2, 50, tanh),  # Dense layer
                    FastDense(50, 2)) # Dense layer without activation

# Build a a neural ode, whenever you want to calculate the derivative call the NN
neural_ode_f(u,p,t) = dudt2(u,p)  

pinit = initial_params(dudt2) # Initial parameters of your ODE come from the intial parameters of your NN

# prob = ODEProblem(neural_ode_f, u0, (0.0f0, 1.5f0), pinit)  # We solve a NN using numerical ODE solvers
prob = ODEProblem(neural_ode_f, u0, tspan, pinit)  # We solve a NN using numerical ODE solvers
sol = solve(prob)

plot(sol)
scatter!(tsteps, ode_data')

# Here we follow the steps for inference of parameters
function loss(p) # This is our loss function
    tmp_prob = remake(prob, p=p)
    tmp_sol = solve(tmp_prob, Tsit5(), saveat = tsteps)
    sum(abs2, Array(tmp_sol) - ode_data)
end

function neuralode_callback(p,l)
    @show l # Shows the lss
    tmp_prob = remake(prob, p=p)
    tmp_sol = solve(tmp_prob, Tsit5(), saveat = tsteps)
    fig = plot(tmp_sol)
    scatter!(fig, tsteps, ode_data')
    display(fig)
    false # A boolean tells you if the optimizaqtion continues or not (false means no, it won't stop)
end

# loss(pinit) # just to print initial error (with 'random parameters')

# Starts with Adam then moves to BFGS
# Search for controlling the adjoints (https://diffeqflux.sciml.ai/dev/ControllingAdjoints/)
@time res = DiffEqFlux.sciml_train(loss, pinit, ADAM(0.05), 
                    maxiters = 100,
                    cb = neuralode_callback)
@time res2 =  DiffEqFlux.sciml_train(loss, res.minimizer, BFGS(initial_stepnorm=0.01), cb = neuralode_callback)


###################### Sensitivity algorithms (2:55:23) #########################
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


###################### Universal ODE (3:04:47) #########################
# The idea here is that you start with 'some part' of your DE and recover the missing terms.
# https://github.com/ChrisRackauckas/universal_differential_equations
# Step 1 lotka
# Step 2 Generate data with the real solution just for the first time steps
# Step 3 Add some noise
# Step 4 Define the UDE
using OrdinaryDiffEq
using ModelingToolkit
using DataDrivenDiffEq
using LinearAlgebra, DiffEqSensitivity, Optim
using DiffEqFlux, Flux
using JLD2


# Simulate dataset (with the solution as an array)
function lotka_volterra!(du, u, p, t)
    r, w = u  # Obtain current values of rabbits and wolves from U
    α, β, γ, δ = p
    du[1] = dr = α*r - β*r*w
    du[2] = dw = γ*r*w - δ*w
end

print("Parameters")
tspan = (0.0, 3.0)
u₀ = [0.44249296, 4.6280594]
p = [1.3, 0.9, 0.8, 1.8] #This will be the TRUE solution

print("Simulating data (true solution)")
prob = ODEProblem(lotka_volterra!, u₀, tspan, p)
sol = solve(prob, saveat=0.1)  # Here you can choose your solver and the accuracy desired options are here: https://diffeq.sciml.ai/dev/basics/common_solver_opts/
# Here is where we get an array from the solution
dataset = Array(sol)
noisydata = dataset + Float32(1e-5)*rand(eltype(dataset), size(dataset))
plot(sol)
# scatter!(sol.t, dataset')  
scatter!(sol.t, noisydata')   

## Simulate dataset (with the solution as an array)
ann = FastChain(FastDense(2,32,tanh), FastDense(32, 32, tanh))
p = initial_params(ann)

# Here we define our UDE. We know some parts but not all of them
function dudt(u, p, t)
    x, y = u
    z = ann(u, p) # 
    [p[1]*x + z[1],
     p[4]*y  + z[2]]
end
    
prob_nn = ODEProblem(dudt, u₀, tspan, p)
sol = concrete_solve(prob_nn, Tsit5(), u₀, p, saveat=0.1)  # Here you can choose your solver and the accuracy desired options are here: https://diffeq.sciml.ai/dev/basics/common_solver_opts/
plot(sol)
scatter!(sol.t, dataset')  # Transposed dataset, just because of the way plots and difeq saves the 

##
function predict(θ)
    Array(concrete_solve(prob_nn, Vern7(), u₀, θ, saveat=0.1, abstol=1e-6, reltol=1e-6))
end

function loss(θ)
    pred = predict(θ)
    sum(abs2, noisydata .- pred), pred
end

loss(p)

const losses = []

function train_callback(θ, l, pred)
    push!(losses, l)
    if length(losses)%50 == 0
        print(losses[end])
    end
    false
end


@time res = DiffEqFlux.sciml_train(loss, p, ADAM(0.01), maxiters = 100, cb = train_callback)
@time res2 = DiffEqFlux.sciml_train(loss, res.minimizer, BFGS(initial_stepnorm=0.01), cb = train_callback)

plot(losses, yaxis = :log, xaxis = :log, xlabel = "Iterations", ylabel = "Loss")

NNsolution = predict(res2.minimizer)
plot(sol.t, noisydata')
scatter!(sol.t, NNsolution')

## ====================== Then we transform our NN back into equations (3:13:30) ========================
using OrdinaryDiffEq
using ModelingToolkit
using DataDrivenDiffEq
using LinearAlgebra, DiffEqSensitivity, Optim
using DiffEqFlux, Flux
using Plots
@variables u[1:2]

# Here you create the symbolic expressions
polys = [] 
for i ∈ 1:5
    push!(polys, u[1]^i)
    push!(polys, u[2]^i)
    for j ∈ 1:5
        if i != j
            push!(polys, (u[1]^i)*(u[2]^j))
            push!(polys, u[2]^i*u[1]^j)
        end
    end
end

# Building a basis of a bunch of polynomials
hh = [cos.(u)...; sin.(u)...; polys...]
basis = Basis(hh,u)

# Optimizer for the SINDY problem
opt = SR3()
# Thresholds which should be used in the search process
λ = exp10.(-6:0.1:2)

f_target(x, w) = iszero(x[1]) ? Inf : norm(w.*x, 2)

X = noisydata
DX = Array(solu)

println("SINDy on full ideal, unavailable data")
Ψ = SINDy(Xₙ[:, :], DX[:, :], basis, λ, opt, g = g, maxiter = 10000) # Fail
println(Ψ)
print_equations(Ψ)

# Test on uode derivative data
println("SINDy on learned, partial, available data")
Ψ = SINDy(Xₙ[:, 2:end], L̂[:, 2:end], basis, λ,  opt, g = g, maxiter = 10000, normalize = true, denoise = true) # Succeed
println(Ψ)
print_equations(Ψ)


# datadriven.sciml.ai/devl/sparse_identification/sindy
# Solving differential equations with neural networks
# # NeuralPDE.jl IMPORTANT
