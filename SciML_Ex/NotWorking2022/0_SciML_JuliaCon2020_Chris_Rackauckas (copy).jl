# Usings
# https://live.juliacon.org/talk/C9FGPP   (Chris workshop)
# https://sciml.ai/
# https://github.com/SciML/SciMLTutorials.jl   (Tutorials)
# https://tutorials.sciml.ai/
# https://github.com/SciML/NeuralPDE.jl    (for PDE)
# https://diffeq.sciml.ai/dev/solvers/ode_solve/   (Here you can see the methods and suggests which one to use)
using DifferentialEquations
using Plots
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
prob = ODEProblem(lotka_volterra!, u₀, tspan, p)
sol = solve(prob)  # Here you can choose your solver and the accuracy desired options are here: https://diffeq.sciml.ai/dev/basics/common_solver_opts/
new_prob = remake(prob, p=[1.2, 1.0, 3.0, 1.0])

# println("Solution values: $sol")
# println("Solution values by index: $(sol[5])")
# println("Solution values by time step: $(sol(.3))")
println("Making some plots.....")
plot(sol)
plot(sol, vars=(1,2)) # Plotting the phase-space (rabbits vs wolves)


## ------------ Example 2 (stochastic ODE)------------ 1:30:53 from video
# Adding a random part into our differential equation
function multiplicative_noise!(du, u, p, t)
    r, w = u
    du[1] = 0.3*r  # This noise depend on the rabbits
    du[2] = 0.3*w  # This noise depend on the wolves
end

# Example 1 Stochastic DE (we add the noise and change the solver to SDE)
prob = SDEProblem(lotka_volterra!, multiplicative_noise!, u₀, tspan, p)
sol = solve(prob)

plot(sol)

## ---------- Example 3 using an ensamble (Ensemble simulations: https://diffeq.sciml.ai/dev/features/ensemble/ )
ensembleprob = EnsembleProblem(prob)
sol = solve(ensembleprob, trajectories =100)
# sol = solve(ensembleprob, SOSRI(), EnsembleThreads(), trajectories =100)
plot(sol)
summ = EnsembleSummary(sol)
plot(summ)


## ------------ Example 4 (delayed ODE 1:39:10 video)------------ https://diffeq.sciml.ai/v2.0/tutorials/dde_example.html
#  You specify a history function h(t) which uses interpolations throught the solution history to form a continuous
# extension of the solver's past.

# In this example we want to add a delay between when the amount of grow on the rabits depends
# on the amount of rabits in a previous time. 
τ = 1.0
function lotka_volterra!(du, u, h, p, t)
    r, w = u  # Obtain current values of rabbits and wolves from U
    del_r = h(p, t - τ; idxs=1)  # idxs means that we want to get the first index
    α, β, γ, δ = p
    du[1] = dr = α*del_r - β*r*w  # Here we use the delayed rabit to decide the grow rate
    du[2] = dw = γ*r*w - δ*w
end

u₀ = [1.0, 1.0]
tspan = (0.0, 10.0)
p = [1.5, 1.0, 3.0, 1.0]
h(p, t) = [1.0, 1.0] 
h(p, t; idxs=1) = 1.0 # Not sure

print("Solving!")
prob = DDEProblem(lotka_volterra!, u₀, h, tspan, p, constant_lag = [τ]) # We define a constant lag tau
sol = solve(prob)
plot(sol)
print("Done!")
prob = ODEProblem(lotka_volterra!, u₀, tspan, p)
sol = solve(prob)
plot(sol)


## ----- Example 5 Calculations control, using callbacks https://diffeq.sciml.ai/dev/features/callback_functions/#Using-Callbacks
# In this example, everytime the numer of wolves is 4 we allow hunting of 1 day (killing '1' amount of wolves)
burn_rabbit_condition(u, t, integrator) = u[2] - 4   # We search for a 0-crossing  (that is why this means u[1] = 0)
burn_rabbit_affect!(integrator) = integrator.u[2] -=1  # Whenever we have '4' wolves, we kill '1'
cb = ContinuousCallback(burn_rabbit_condition, burn_rabbit_affect!)
prob = DDEProblem(lotka_volterra!, u₀, h, tspan, p, constant_lag = [τ]) # We define a constant lag tau
sol = solve(prob, callback = cb)
plot(sol)

## ----- Example 6 algebraic equations https://diffeq.sciml.ai/dev/tutorials/advanced_ode_example/
# How to incorporate restrictions into your model, like a concervation law


## 
# Simulate dataset (with the solution)

dataset = Array(sol)
scatter!(sol.t, dataset')

# Changing the paramteres
tmp_problem = remake(prob, u\_0, p = [1.2, 0.8, 2.5, 0.8])

tmp_solve = solve(tmp_prob)
plot(tmp_solve)
scatter!(sol.t, dataset')

function loss(p)
    tmp_prob = remake(prob, p=p)
    tmp_sol = solve(tmp_prob, saveat=0.1, verbose=false)
    sum(abs2, Array(tmp_sol) - dataset), tmp_sol
end

using DiffEqFlux, Optim

# Initial guess
pinit = [1.2, 0.8, 2.5, 0.8]
res = DiffEqFlux.sciml_train(loss, BFGS())

# What if found
res.minimizer

function plot_callback(p, l, tmp_sol)
    @show l
    tmp_prob = remake(prob, p=p)
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

prob = ODEProblem(neural_ode_f, u0, (0.0f0, 1.5f0), pinit)
sol = solve(prob)

plot(sol)
scatter!(tsteps, ode_data')

function loss(p)
    tmp_prob = remake(prob, p=p)
    tmp_sol = solve(tmp_prob, Tsit5(), saveat = tsteps)
    sum(abs2, Array(tmp_sol) - ode_data)
end

function neuralode_callback(p,l)
    @show l
    tmp_prob = remake(prob, p=p)
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
    tmp_prob = remake(prob, p=p)
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
