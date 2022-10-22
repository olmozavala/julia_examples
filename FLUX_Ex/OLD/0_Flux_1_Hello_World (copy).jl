# Examples at: https://github.com/FluxML/model-zoo/
using Flux
using Plots
using Printf
gr() # Using GR to plot

## --------- Make synthetic data
n = 100
x = rand(1,n)*2*π # READ: X and Y MUST be matrices and not vector to work  on the training phase
# y = sin.(x) + rand(length(x))*.1  # Simulated data
# y(x) = 3x + 2 + rand(1)[1]
y(x) = 1.5343x + 2 + rand(1)[1]
x_train = x
y_train = y.(x)

scatter(x', y.(x)', label="Data")

## ------------ Making your model
println("======== Simplest Dense model =================")
# https://fluxml.ai/Flux.jl/v0.4/models/layers.html#Flux.Dense
model = Dense(1 => 1) # One input and one output (by default σ is the identity)
print("Weights: $(model.weight) and bias: $(model.bias)")
plot!(x', model(x)', label="model", c=:red, aspect_ratio=:equal, xlim=[-1,8])

## ------------ Loss and train
loss(a,b) = Flux.Losses.mse(model(a), b)
print("Loss value $(loss(x_train, y_train))") # Current value of loss function (X and Y must be matrices like (1,n) )

##
opt = Descent()
params = Flux.params(model)  # Indicate which parameters we want to consider
mydata = [(x_train, y_train)]
for epoch in 1:200
  Flux.train!(loss, params, mydata, opt)  # Loss function, parameters, data and optimizer
  println("Loss value $(loss(x,y.(x)))") # Current value of loss function (X and Y must be matrices like (1,n) )
end
scatter(x', y.(x)', label="Data",  aspect_ratio=:equal)
plot!(x', model(x)', label="model", c=:red, xlim=[-1,8])

##
# Define loss function 
function loss(x, y)
  ŷ = linmodel(x)
  sum((y .- ŷ).^2)
end

iter = 30
α = 0.02
for i = 1:iter
# Compute gradient of W and b with respect to the loss
  gs = gradient(() -> loss(x,y), params(W, b))
  W .-= α .* gs[W]
  b .-= α .* gs[b]
  println("After gds loss: $(loss(x,y)) W: $W and b: $b")
  display(plot!(x, linmodel(x), label="Pred$i"))
  sleep(.1)
end
scatter(x, y, label="Data")
plot!(x, linmodel(x), label="Last")
println("W=$W b=$b")

println("======== Composed models =================")

## -------------------
using Flux

actual(x) = 4x + 2
n = 10

x_train = rand(1,n).*2*π
y_train = actual.(x_train)

scatter(x_train', y_train')

predict = Dense(1 => 1)
loss(x, y) = Flux.Losses.mse(predict(x), y);
opt = Descent()
parameters = Flux.params(predict)
data = [(x_train, y_train)]
for epoch in 1:10
  Flux.train!(loss, parameters, data, opt)
  println(loss(x_train, y_train))
end