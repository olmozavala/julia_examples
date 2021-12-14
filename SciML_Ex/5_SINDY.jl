## Usings
# https://live.juliacon.org/talk/C9FGPP   (Chris workshop)
# https://sciml.ai/
# https://github.com/SciML/SciMLTutorials.jl   (Tutorials)
# https://tutorials.sciml.ai/
# https://github.com/SciML/NeuralPDE.jl    (for PDE)
# https://diffeq.sciml.ai/dev/solvers/ode_solve/   (Here you can see the methods and suggests which one to use)
using DifferentialEquations
using Plots
using OrdinaryDiffEq
using DataDrivenDiffEq
using LinearAlgebra, DiffEqSensitivity, Optim
using DiffEqFlux, Flux

## #################### Universal ODE (3:04:47) #########################
# The idea here is that you start with 'some part' of your DE and recover the missing terms.
# https://github.com/ChrisRackauckas/universal_differential_equations

# --------------------- Simulate dataset (with the solution as an array) and plot the training data ------------
function lotka_volterra!(du, u, p, t)
    r, w = u  # Obtain current values of rabbits and wolves from U
    α, β, γ, δ = p
    du[1] = dr = α*r - β*r*w
    du[2] = dw = γ*r*w - δ*w
end

print("Parameters")
tspan = (0.0, 3.0)
u0 = Float32[0.44249296,4.6280594]
p_ = Float32[1.3, 0.9, 0.8, 1.8]

print("Simulating data (true solution)")
prob = ODEProblem(lotka_volterra!, u₀, tspan, p_orig)
sol_orig = solve(prob, saveat=0.1)  # Here you can choose your solver and the accuracy desired options are here: https://diffeq.sciml.ai/dev/basics/common_solver_opts/
# Here is where we get an array from the solution
dataset = Array(sol_orig)
noisydata = dataset + Float32(1e-5)*rand(eltype(dataset), size(dataset))
plot(sol_orig)
# scatter!(sol.t, dataset')  
scatter!(sol_orig.t, noisydata')   

## ------------ Make your NN and create an ODE where two terms are represented as a NN
ann =  FastChain(FastDense(2, 32, tanh),FastDense(32, 32, tanh), FastDense(32, 2))

# Here we define our UDE. We know some parts but not all of them
function dudt(u, p, t)
    x, y = u
    z = ann(u, p) # Here we define z as aour NN
    [p_orig[1]*x + z[1],
     z[2] - p_orig[4]*y ]
end
    
p = initial_params(ann)
prob_nn = ODEProblem(dudt, u₀, tspan, p)
sol = concrete_solve(prob_nn, Tsit5(), u₀, p, saveat=0.1)  # Here you can choose your solver and the accuracy desired options are here: https://diffeq.sciml.ai/dev/basics/common_solver_opts/
plot(sol)
scatter!(sol.t, dataset', title="Solution with initial parameters")  # Transposed dataset, just because of the way plots and difeq saves the 

## === Now train the NN parameters inside the ODE

# Just a function to solve with a Numerical Method (RK4 for example)
function predict(θ)
    Array(concrete_solve(prob_nn, Vern7(), u₀, θ, saveat=0.1, abstol=1e-6, reltol=1e-6))
end
# Our loss function is the predicted nuerical solution vs our data
function loss(θ)
    pred = predict(θ)
    sum(abs2, noisydata .- pred), pred
end

const losses = []
# Just a callback to print some losses
function train_callback(θ, l, pred)
    push!(losses, l)
    if length(losses)%50 == 0
        println(losses[end])
    end
    false
end

# Here we optimize the parameters of the NN part
println("Running the training...")
println("ADAM ...")
@time res = DiffEqFlux.sciml_train(loss, p, ADAM(0.01), maxiters = 100, cb = train_callback)
println("BFGS ...")
@time res2 = DiffEqFlux.sciml_train(loss, res.minimizer, BFGS(initial_stepnorm=0.01), cb = train_callback)
println("Done! .... making some  final plots")

plot(losses, yaxis = :log, xaxis = :log, xlabel = "Iterations", ylabel = "Loss")

# NNsolution = predict(res2.minimizer) # Here we are solving for the original tspan 
NNsolution = Array(concrete_solve(prob_nn, Vern7(), u₀, res2.minimizer, saveat=0.1, abstol=1e-6, reltol=1e-6))
plot(sol.t, dataset', title="Data and solved ODE with NN inside of it", label= ["Original u₁" "Original u₂"])
scatter!(sol.t, NNsolution', label= ["NN u₁" "NN u₂"])

##  Show the values of the derivative of both variables u₀ and u₁ for the original and computed solutions

# https://docs.julialang.org/en/v1/base/base/#Base.Val
DX = Array(sol_orig(sol.t, Val{1})) # Obtaining the numerical derivative of the solution 

prob_nn2 = ODEProblem(dudt, u₀, tspan, res2.minimizer) # Solution with current obtained parameters
sol_nn2 = solve(prob_nn2, Tsit5())  # Solve
DX_nn2 = Array(sol_nn2(sol.t, Val{1})) # Obtaining the numerical derivative of the new solution

# Plotting the learned derivatives
plot(DX')
plot!(DX_nn2')

# ==========================================================================================
# 
## ====================== Then we transform our NN back into equations (3:13:30) ========================
@variables u[1:2]

# Here you create the symbolic expressions
using ModelingToolkit

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
println("The obtained basis is: \n", basis)
## Now we test a couple versions of SINDY

# Optimizer for the SINDY problem
opt = SR3()
# Thresholds which should be used in the search process
λ = exp10.(-6:0.1:2)

# g(x) = x[1] < 1 ? Inf : norm(x, 2)
g(x, w) = iszero(x[1]) ? Inf : norm(w.*x, 2)

println("SINDy on full ideal, unavailable data")
Ψ = SINDy(noisydata[:,:], DX[:, :], basis, λ, opt,  maxiter = 10000, g = g) # Fail
println(Ψ)
print_equations(Ψ)

## ======== Now we try only using data from the missing terms (3:17)
    # du[1] = dr = α*r - β*r*w
    # du[2] = dw = γ*r*w - δ*w
L = [-p_orig[2]*(X[1,:].*X[2,:])';p_orig[3]*(X[1,:].*X[2,:])']  # These are the terms that the nn wants to match (without noise)
# Neural network guess
L̂ = ann(noisydata,res2.minimizer)

scatter(L', label=["True u₁", "True u₂" ])
plot!(L̂', label=["NN u₁", "NN u₂" ])

# Plot the error
scatter(abs.(L-L̂)', yaxis = :log, title="Error")
println("SINDy on learned, partial, available data")
# Ψ = SINDy(noisydata[:, 2:end], L̂[:, 2:end], basis, λ,  opt, g = g, maxiter = 10000, normalize = true, denoise = true) # Succeed
Ψ = SINDy(noisydata[:, 2:end], L̂[:, 2:end], basis, λ,   g = g,  normalize = true, denoise = true) # Succeed
println(Ψ)
print_equations(Ψ)


# # datadriven.sciml.ai/devl/sparse_identification/sindy
# # Solving differential equations with neural networks
# # # NeuralPDE.jl IMPORTANT
