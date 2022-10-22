using DifferentialEquations
using Plots
using OrdinaryDiffEq
using ModelingToolkit
using DataDrivenDiffEq
using LinearAlgebra, Optim
using DiffEqFlux, Flux

## Simulate dataset (with the solution as an array)
function lotka_volterra!(du, u, p, t)
    r, w = u  # Obtain current values of rabbits and wolves from U
    α, β, γ, δ = p
    du[1] = dr = α*r - β*r*w
    du[2] = dw = γ*r*w - δ*w
end

println("Parameters")
tspan = (0.0, 5.0)
u₀ = [1., 1.]
p = [1.3, 0.9, 0.8, 1.8] #This will be the TRUE solution
original_p = p

println("Simulating data (true solution)....")
prob = ODEProblem(lotka_volterra!, u₀, tspan, original_p)
sol = solve(prob, saveat=0.1)  # Here you can choose your solver and the accuracy desired options are here: https://diffeq.sciml.ai/dev/basics/common_solver_opts/
dataset = Array(sol)
# noisydata = dataset + Float32(1e-1)*rand(eltype(dataset), size(dataset))
noisydata = dataset 
plot(sol)
scatter!(sol.t, noisydata', title="Example Data")

## Make your Partial NODE
ann = FastChain(FastDense(2,32,tanh), FastDense(32, 32, tanh), FastDense(32, 2))
pnn = initial_params(ann)

# Here we define our UDE. We know some parts but not all of them
function dudt(u, p, t)
    x, y = u
    z = ann(u, p) # Here we define z as aour NN
    [original_p[1]*x + z[1],
     original_p[4]*y  + z[2]]
end
    
prob_nn = ODEProblem(dudt, u₀, tspan, pnn)
sol = concrete_solve(prob_nn, Tsit5(), u₀, pnn, saveat=0.1)  # Here you can choose your solver and the accuracy desired options are here: https://diffeq.sciml.ai/dev/basics/common_solver_opts/
plot(sol)
scatter!(sol.t, dataset', title="Solution with initial parameters")  # Transposed dataset, just because of the way plots and difeq saves the 

## === Now optimizing NN parameters inside the ODE
# Just a function to solve with a Numerical Method (RK4 for example)
function predict(θ)
    Array(concrete_solve(prob_nn, Vern7(), u₀, θ, tspan=tspan, saveat=0.1, abstol=1e-6, reltol=1e-6))
end
# Our loss function is the predicted nuerical solution vs our data
function loss(θ)
    pred = predict(θ)
    sum(abs2, noisydata .- pred), pred
end

## Just before training...
const losses = []
# Just a callback to print some losses
function train_callback(θ, l, pred)
    push!(losses, l)
    if length(losses)%50 == 0
        println(losses[end])
    end
    false
end

## Here we optimize the parameters of the NN part
@time res = DiffEqFlux.sciml_train(loss, pnn, ADAM(0.01), maxiters = 100, cb = train_callback)
println("Final error Adam of: $(res.minimum)")
@time res2 = DiffEqFlux.sciml_train(loss, res.minimizer, BFGS(initial_stepnorm=0.01), cb = train_callback)
println("Final error of: $(res2.minimum)")
plot(losses, yaxis = :log, xaxis = :log, xlabel = "Iterations", ylabel = "Loss")

## Plotting the results
NNsolution = predict(res2.minimizer)
plot(sol.t, dataset', title="Data and solved ODE with NN inside of it", label= ["Original u₁" "Original u₂"])
scatter!(sol.t, NNsolution', label= ["NN u₁" "NN u₂"])

## Plot the derivatives found by the NN
prob_nn2 = ODEProblem(dudt, u₀, tspan, res2.minimizer)
sol_nn = solve(prob_nn2, saveat=0.1)

X = noisydata
DX = Array(sol(sol.t, Val{1}))   # Gets DX
DX_nn = Array(sol_nn(sol.t, Val{1}))

plot(DX')
plot!(DX_nn')

## Now lets predict outside the training domain
test_tspan = (0.0, 4.0)

test_prob = ODEProblem(lotka_volterra!, u₀, test_tspan, original_p)
test_prob_nn = ODEProblem(dudt, u₀, test_tspan, res2.minimizer)

test_sol = solve(test_prob, saveat=0.1)  # Here you can choose your solver and the accuracy desired options are here: https://diffeq.sciml.ai/dev/basics/common_solver_opts/
test_sol_nn = solve(test_prob_nn, saveat=0.1)
plot(test_sol, label= ["Original u₁" "Original u₂"])
scatter!(test_sol.t, test_sol_nn', label= ["NN u₁" "NN u₂"])
vline!([3.0])


## ====================== Then we transform our NN back into equations (3:13:30) ========================
# using OrdinaryDiffEq
# using ModelingToolkit
# using DataDrivenDiffEq
# using LinearAlgebra, DiffEqSensitivity, Optim
# using DiffEqFlux, Flux
# using Plots
# @variables u[1:2]

# # Here you create the symbolic expressions
# polys = [] 
# for i ∈ 1:5
#     push!(polys, u[1]^i)
#     push!(polys, u[2]^i)
#     for j ∈ 1:5
#         if i != j
#             push!(polys, (u[1]^i)*(u[2]^j))
#             push!(polys, u[2]^i*u[1]^j)
#         end
#     end
# end

# # Building a basis of a bunch of polynomials
# hh = [cos.(u)...; sin.(u)...; polys...]
# basis = Basis(hh,u)

# # Optimizer for the SINDY problem
# opt = SR3()
# # Thresholds which should be used in the search process
# λ = exp10.(-6:0.1:2)

# f_target(x, w) = iszero(x[1]) ? Inf : norm(w.*x, 2)

# X = noisydata
# DX = Array(solu)

# println("SINDy on full ideal, unavailable data")
# Ψ = SINDy(Xₙ[:, :], DX[:, :], basis, λ, opt, g = g, maxiter = 10000) # Fail
# println(Ψ)
# print_equations(Ψ)

# # Test on uode derivative data
# println("SINDy on learned, partial, available data")
# Ψ = SINDy(Xₙ[:, 2:end], L̂[:, 2:end], basis, λ,  opt, g = g, maxiter = 10000, normalize = true, denoise = true) # Succeed
# println(Ψ)
# print_equations(Ψ)


# # datadriven.sciml.ai/devl/sparse_identification/sindy
# # Solving differential equations with neural networks
# # # NeuralPDE.jl IMPORTANT
