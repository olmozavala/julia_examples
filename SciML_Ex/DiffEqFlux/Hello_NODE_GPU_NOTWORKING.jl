using Flux, DiffEqFlux, Optimization, OptimizationFlux, Zygote, 
      OrdinaryDiffEq, Plots, CUDA, SciMLSensitivity, Random, ComponentArrays
CUDA.allowscalar(false) # Makes sure no slow operations are occuring

#rng for Lux.setup
rng = Random.default_rng()
# Generate Data
datasize = 30
tspan = (0.0f0, 1.5f0)
tsteps = range(tspan[1], tspan[2], length = datasize)
u0 = Float32[0.44249296,4.6280594]
p_ = Float32[1.3, 0.9, 0.8, 1.8]


function lotka!(du, u, p, t)
    α, β, γ, δ = p
    du[1] = α*u[1] - β*u[2]*u[1]
    du[2] = γ*u[1]*u[2]  - δ*u[2]
end
p = Float32[1.3, 0.9, 0.8, 1.8]

prob_trueode = ODEProblem(lotka!, u0, tspan, p)

# Make the data into a GPU-based array if the user has a GPU
ode_data = gpu(solve(prob_trueode, Tsit5(), saveat = tsteps))
plot(ode_data)

## -----------------
# dudt2 = Chain(x -> x.^3, Dense(2, 50, tanh), Dense(50, 2)) |> gpu
dudt2 = Chain(Dense(2, 50, tanh), Dense(50, 2)) |> gpu
u0 = Float32[2.0; 0.0] |> gpu
prob_neuralode = NeuralODE(dudt2, tspan, Tsit5(), saveat = tsteps)

function predict_neuralode(p)
  gpu(prob_neuralode(u0,p))
end
function loss_neuralode(p)
    pred = predict_neuralode(p)
    loss = sum(abs2, ode_data .- pred)
    return loss, pred
end
# Callback function to observe training
list_plots = []
iter = 0
callback = function (p, l, pred; doplot = false)
  global list_plots, iter
  if iter == 0
    list_plots = []
  end
  iter += 1
  display(l)
  # plot current prediction against data
  plt = scatter(tsteps, Array(ode_data[1,:]), label = "data")
  scatter!(plt, tsteps, Array(pred[1,:]), label = "prediction")
  push!(list_plots, plt)
  if doplot
    display(plot(plt))
  end
  return false
end

adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x,p)->loss_neuralode(x), adtype)
optprob = Optimization.OptimizationProblem(optf, prob_neuralode.p)
result_neuralode = Optimization.solve(optprob,ADAM(0.05),callback = callback,maxiters = 300)