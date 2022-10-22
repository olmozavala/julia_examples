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
# 🐰  🔥 

## ------------------ First example classic lotka_volterra ODE
# Time 1:17 of Chris talk
function lotka_volterra!(du, u, p, t)
    r, w = u  # Obtain current values of rabbits and wolves from U
    α, β, γ, δ = p
    du[1] = dr = α*r - β*r*w
    du[2] = dw = γ*r*w - δ*w
end

u₀ = [1.0, 1.0]
tspan = (0.0, 20.0)
p = [1.2, 0.6, 0.3, 0.8]

print("Solving the problem....")
orig_prob = ODEProblem(lotka_volterra!, u₀, tspan, p)
sol = solve(orig_prob)  # Here you can choose your solver and the accuracy desired options are here: https://diffeq.sciml.ai/dev/basics/common_solver_opts/
# new_prob = remake(orig_prob, p=[1.2, 0.6, 0.3, 0.8])

# println("Solution values: $sol")
# println("Solution values by index: $(sol[5])")
# println("Solution values by time step: $(sol(.3))")
println("Making some plots.....")
plot(sol)
# 
## 
# Simulate dataset (with the solution)

dataset = Array(sol)
scatter!(sol.t, dataset')

# Changing the paramteres
# tmp_solve = solve(new_prob)
# plot(tmp_solve)
# scatter!(sol.t, dataset')

function loss(p)
    tmp_prob = remake(orig_prob, p=p)
    tmp_sol = solve(tmp_prob, saveat=0.1, verbose=false)
    sum(abs2, Array(tmp_sol) - dataset), tmp_sol
end

# Initial guess
pinit = [1.2, 0.8, 2.5, 0.8]
res = DiffEqFlux.sciml_train(loss, BFGS())

## What if found
res.minimizer

function plot_callback(p, l, tmp_sol)
    @show l
    tmp_prob = remake(orig_prob, p=p)
    tmp_sol = solve(tmp_prob, saveat=0.1, verbose=false)
    fig = plot(tmp_sol)
    scatter(sol.t, dataset')
    display(fig)
    false
end


using flux
res = DiffEqFlux.sciml_train(loss, pinit, ADAM(0.01), cb=plot_callback, maxiters=1000)
res.minimizer

using BlackBoxOptim

res = DiffEqFlux.sciml_train(loss, pinit, DiffEqFlux.BBO(), lower_bownds=0.5ones(4), upper_bounds=4.0ones(4))

# Deine the same model

Turing.setadbackend(:forwarddiff)

@model funciton fitlv(data)
    σ ~ InverseGamma(2, 3)
    α ~ truncated(Normal(1.3, 0.5), 0.5, 2.5)
    β ~ truncated(Normal(1.2, 0.5), 0, 2.5)
    α ~ truncated(Normal(1.3, 0.5), 0.5, 2.5)
    α ~ truncated(Normal(1.3, 0.5), 0.5, 2.5)

    α ~ truncated(Normal(1.3, 0.5), 0.5, 2.5)


using DiffEqFlux, OrdinaryDiffEq, Flux, Optim, Plots

u0 = Float32[2.0; 0.0]
datasize = 30
tspan = (0.0f0, 1.5f0)
tsteps = range(tspan[1], tspan[2], length = datasize)

function trueODEfunc(du, u, p, t)
    true_A = [-0.1 2.0; -2.0 -0.1]
    du .= ((u.^3)'true_a)'
end

prob_trueode = ODEProblem(trueODEfunc, u0, tspan)
ode_data = Array(solve(prob_trueode, Tsit5(), saveat = tsteps))
    
dudt2 = FastChain((x, p) -> x.^3,
                    FastDense(2, 50, tanh), 
                    FastDense(50, 2))

neural_ode_f(u,p,t) = dudt2(u,p)

pinit = initial_params(dudt2)

orig_prob = ODEProblem(neural_ode_f, u0, (0.0f0, 1.5f0), pinit)
sol = solve(orig_prob)

plot(sol)
scatter!(tsteps, ode_data')

function loss(p)
    tmp_prob = remake(orig_prob, p=p)
    tmp_sol = solve(tmp_prob, Tsit5(), saveat = tsteps)
    sum(abs2, Array(tmp_sol) - ode_data)
end

function neuralode_callback(p,l)
    @show l
    tmp_prob = remake(orig_prob, p=p)
    tmp_sol = solve(tmp_prob, Tsit5(), saveat = tsteps)
    fig = plot(tmp_sol)
    scatter!(fig, tsteps, ode_data')
    display(fig)
    false
end

loss(pinit)

# Starts with Adam then moves to BFGS
@time res = DiffEqFlux.sciml_train(loss, pinit, ADAM(0.05), 
                    maxiters = 100,
                    cb = neuralode_callback)

@time res2 =  DiffEqFlux.sciml_train(loss, res, BFGS(initial_stepnorm=0.01) cb = neuralode_callback)

# Search for controlling the adjoints

function loss(p)
    tmp_prob = remake(orig_prob, p=p)
    tmp_sol = solve(tmp_prob, Tsit5(), saveat = tsteps)
    sum(abs2, Array(tmp_sol) - ode_data)
end


# Step 1 lotka

# Step 2 Generate data with the real solution just for the first time steps

# Step 3 Add some noise

# Step 4 Define the UDE



ann = FastChain(FasteDense(2,e2,tanh), FastDense(32, 32, tanh))

function dudt_(u, p, t)
    x,y = u
    z = ann(u, p) # 
    [p_p[1]*x + z[1]],
    _p_[4]*y 
    


@variables u[1:2]
polys = 

datadriven.sciml.ai/devl/sparse_identification/sindy



# olving differential equations with neural networks


# NeuralPDE.jl IMPORTANT
