# Examples at: https://github.com/FluxML/model-zoo/
using Flux
using Plots
using Printf
using Statistics
gr() # Using GR to plot

## --------- Make synthetic data
n = 100
x = rand(1,n)*2π # READ: X and Y MUST be matrices and not vector to work  on the training phase
# y(x) = 5x + 2 + rand(1)[1]*3
y(x) = sin(x) + rand(1)[1]*.5
# Normalize data
x_train = (x .- mean(x)) ./ std(x)
y_train = y.(x)
y_train = (y_train .- mean(y_train))./ std(y_train)

## ------------ Making your model
println("======== Simplest Dense model =================")
# https://fluxml.ai/Flux.jl/v0.4/models/layers.html#Flux.Dense

model = Dense(1 => 1) # One input and one output (by default σ is the identity)
# model = Chain(Dense(1 => 1)) # One input and one output (by default σ is the identity)
act = relu
model = Chain(Dense(1 => 5, act), 
              Dense(5 => 5, act), 
              Dense(5 => 1)) 

print("Weights: $(model.weight) and bias: $(model.bias)")
scatter(x_train', y_train', label="Data")
plot!(x_train', model(x_train)', label="model", c=:red)

loss(a,b) = Flux.Losses.mse(model(a), b) # Loss function 

opt = Adam(0.01)

params = Flux.params(model)  # Indicate which parameters we want to consider
mydata = [(x_train, y_train)]
x_sort = sort(x_train, dims = 2)
@gif for epoch in 1:500
  Flux.train!(loss, params, mydata, opt)  # Loss function, parameters, data and optimizer
  if epoch % 10 == 0
    scatter(x_train', y_train', label="Data")
    p = plot!(x_sort', model(x_sort)', label="model", c=:red, title="Epoch $epoch")
    println("Loss value $(loss(x,y_train))") # Current value of loss function (X and Y must be matrices like (1,n) )
  end
end

##
scatter(x_train', y_train', label="Data")
plot!(x_sort', model(x_sort)', label="model", c=:red, title="Final")