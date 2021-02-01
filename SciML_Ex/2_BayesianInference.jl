# Usings
# https://live.juliacon.org/talk/C9FGPP   (Chris workshop)
# https://sciml.ai/
# https://github.com/SciML/SciMLTutorials.jl   (Tutorials)
# https://tutorials.sciml.ai/
# https://github.com/SciML/NeuralPDE.jl    (for PDE)
# https://diffeq.sciml.ai/dev/solvers/ode_solve/   (Here you can see the methods and suggests which one to use)
using DiffEqFlux, OrdinaryDiffEq, Flux, Optim, Plots, Turing, MCMCChains, StatsPlots

## ----------------------------------------------------------------------------------------
## --------------Bayessian inference (predict parameters with uncertainty) (2:23 from workshop) ------------------------------------------------
# Website: https://turing.ml/dev/tutorials/10-bayesiandiffeq/     
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


