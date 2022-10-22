# Examples at: https://github.com/FluxML/model-zoo/
using Flux
using Plots
using Printf
using Statistics
using CUDA
CUDA.functional()
gr() # Using GR to plot

## --------- Make synthetic data
n = 100
x = rand(1,n)*2π # READ: X and Y MUST be matrices and not vector to work  on the training phase
x = sort(x, dims = 2)
# y(x) = 5x + 2 + rand(1)[1]*3
y(x) = sin(x) + rand(1)[1]*.5
# Normalize data
x_train_cpu = (x .- mean(x)) ./ std(x)
y_train_cpu = y.(x)
y_train_cpu = (y_train_cpu .- mean(y_train_cpu))./ std(y_train_cpu)
x_train = cu(x_train_cpu)
y_train = cu(y_train_cpu)
print(typeof(x_train))
scatter(x_train' |> cpu, y_train' |> cpu, label="Data")

## ------------ Making your model
println("======== Simplest Dense model =================")
# ========= Send model weights to GPU with fmap
# model = Chain(Dense(1 => 1)) |> gpu # One input and one output (by default σ is the identity)
act = relu
model = Chain(Dense(1 => 5, act), 
              Dense(5 => 5, act), 
              Dense(5 => 1))  |> gpu

scatter(x_train' |> cpu, y_train' |> cpu, label="Data")
plot!(x_train' |> cpu, model(x_train)' |> cpu, label="model", c=:red)
print("Type $(typeof(model[1].weight))")


loss(a,b) = Flux.Losses.mse(model(a), b) # Loss function 

opt = Adam(0.01)

params = Flux.params(model)  # Indicate which parameters we want to consider
mydata = [(x_train, y_train)]
println("Loss value $(loss(x_train,y_train))") # Current value of loss function (X and Y must be matrices like (1,n) )

println("Running on the GPU...")
@gif for epoch in 1:500
  Flux.train!(loss, params, mydata, opt)  # Loss function, parameters, data and optimizer
  if epoch % 10 == 0
    scatter(x_train' |> cpu, y_train' |> cpu, label="Data")
    p = plot!(x_train' |> cpu, model(x_train)' |> cpu, label="model", c=:red, title="Epoch $epoch")
    println("Loss value $(loss(x_train,y_train))") # Current value of loss function (X and Y must be matrices like (1,n) )
  end
end

##
scatter(x_train', y_train', label="Data")
plot!(x_train', model(x_train)', label="model", c=:red, title="Final")