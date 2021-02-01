# Examples at: https://github.com/FluxML/model-zoo/
using Flux

## Basics gradients
println("======== Basics: Gradients ========")
f(x) = 3*x^2 + 2x + 1
df(x) = gradient(f, x)  # df/dx = 6x + 2

println("Gradient of f $(df(2))")

# Multiple parameters
f(x,y) = sum((x .- y).^2)

println("Gradient of f(x,y) $(gradient(f, [2 ,1], [2, 0]))")  # 4 parameters

# Without passing the parameters
x = [2 ,1]
y = [2 ,0]

gs = gradient(params(x,y)) do 
    sum((x .- y).^2)
end
# df/dx and df/dy
println("df/dx $(gs[x])")
println("df/dy $(gs[y])")

## Simple models
println("======== Simplest model =================")
using Plots
gr() # Using GR to plot

# Random Initial weigths
W = rand(1)   
b = rand(1)

predict(x) = W.*x .+ b

# Define loss function
function loss(x, y)
  ŷ = predict(x)
  sum((y .- ŷ).^2)
end

# Random 'real' data
x = rand(10)
y = 3 .* x .+ 8 + rand(length(x))*.1  # Simulated data

scatter(x, y, label="Data")
plot!(x, predict(x), label="Pred0")
println("Initial loss: $(loss(x, y))") # ~ 3
iter = 10
α = 0.01
for i = 1:iter
# Compute gradient of W and b with respect to the loss
  gs = gradient(() -> loss(x,y), params(W, b))
  W .-= α .* gs[W]
  b .-= α .* gs[b]
  println("After gds loss: $(loss(x,y)) W: $W and b: $b")
  plot!(x, predict(x), label="Pred$i")
end
plot!(x, predict(x), label="Last")


println("======== Composed models =================")