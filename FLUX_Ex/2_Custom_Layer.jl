# Examples at: https://github.com/FluxML/model-zoo/
using Flux
using Plots
using Printf
using Statistics
gr() # Using GR to plot


## --------- Make synthetic data
n = 100
x = rand(1,n)*2π # READ: X and Y MUST be matrices and not vector to work  on the training phase
x = sort(x, dims = 2)
y(x) = 4*sin(x) + 2 + rand(1)[1]*3
# Normalize data
x_train = (x .- mean(x)) ./ std(x)
y_train = y.(x)
y_train = (y_train .- mean(y_train))./ std(y_train)

scatter(x_train', y_train', label="Data")


## ========== Custom Layer =========
struct cSine # Define the parameters in a struct
  a
  b
  c
end
# Allow for tracking params for backprop
Flux.@functor cSine
# Define how are you building the layers
cSine(a::Integer, b::Integer, c::Integer) = cSine(randn(a), randn(b), rand(c))
(m::cSine)(x) = m.a .* sin.(x.*m.b) .+ m.c # Overload call to use it as function
model = cSine(1,1,1)   # One input and output
model(x_train)
# If you call cSine(3,5) is initializing the parameters, if you call cSine(5.313) is calling the model
## ========== Done Custom Layer =========
model = cSine(1,1,1)   # One input and output

scatter(x_train', y_train', label="Data")
plot!(x_train', model(x_train)', label="model", c=:red)

loss(x,y) = Flux.Losses.mse(model(x), y) # Loss function 

## ---------- Training ----------
opt = Adam(0.1)

params = Flux.params(model)  # Indicate which parameters we want to consider
mydata = [(x_train, y_train)]
x_sort = sort(x_train, dims = 2)

@gif for epoch in 1:200
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