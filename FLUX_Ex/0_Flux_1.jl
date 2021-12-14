# Examples at: https://github.com/FluxML/model-zoo/
using Flux

## Basics gradients
println("======== Basics: Gradients ========")
f(x) = 3*x^2 + 2x + 1
df(x) = gradient(f, x)  # df/dx = 6x + 2

println("Gradient of f $(df(2))")

## Multiple parameters
println("======== Basics: Multiple parameters ========")
f(x,y) = sum((x .- y).^2) # Gradient 4 params = (∂f/∂x₁, ∂f/∂x₂, ∂f/∂y₁, ∂f/∂y₂)
println("Gradient of f(x,y) $(gradient(f, [2 ,1], [2, 0]))")  # 4 parameters

## Without passing the parameters
println("======== Basics: Multiple parameters ========")
x = [1,2] # [2 ,1] 

gs = gradient(params(x,y)) do 
    sum((x .- y).^2)
end
# df/dx and df/dy
println("(df/dx, df/dy) = $(gs[x])")  # df/dx sum(2*(x-y)*1)  df/dy (for y  = 2) = sum(2*(x-y)*-1)????

## Simple models
println("======== Simplest Manual model =================")
using Plots
gr() # Using GR to plot

# Random Initial weigths
W = rand(1)   
b = rand(1)

linmodel(x) = W.*x .+ b

# Define loss function 
function loss(x, y)
  ŷ = linmodel(x)
  sum((y .- ŷ).^2)
end

# Random 'real' data
x = rand(10)
y = 3 .* x .+ 8 + rand(length(x))*.1  # Simulated data

scatter(x, y, label="Data")
plot!(x, linmodel(x), label="Pred0")
println("Initial loss: $(loss(x, y))") # ~ 3
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