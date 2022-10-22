
using Flux # MXNet or TensorFlow
using Flux, Statistics
using Flux.Data: DataLoader
using Flux:  @epochs
using Parameters: @with_kw
using CUDA
using Images
using LinearAlgebra
using Plots


x = collect(-π:.1:π)
y = x .+ rand(length(x))*.1
x = reshape(x, 1, 63)
y = reshape(y, 1, 63)
scatter(x', y')

## Define model (simple two layer dense)
mymodel = Chain(
 	    Dense(1, 10, relu),
        Dense(10, 1)
        )


trainmode!(model,true)
loss(x, y) = Flux.mae(mymodel(x), y)
opt = ADAM(0.5)
ps = params(model)

mycb = () -> @show(loss(x,y))

train_loader = Flux.Data.DataLoader((x,y),batchsize=10)
Flux.train!(loss, ps, train_loader, opt, mycb)

scatter(x, y, label="data")
plot!(x, mymodel(x')', label="model")
## 
model = Flux.Chain(Dense(1,1,σ))
a = rand([0,1,2,3,4,5,6,7,8,9], 1,50000)
b = rand([0,1,2,3,4,5,6,7,8,9], 1,50000)
accuracy(x, y, model) = mean((model(x)) .== y)
train_loader = Flux.Data.DataLoader((a,b),batchsize=10)
trainmode!(model,true)
loss(x, y) = Flux.mae(model(x), y)
opt = ADAM(0.5)
ps = params(model)
# Flux.train!(loss, ps, train_loader , opt)
Flux.train!(loss, ps, (a,b), opt)