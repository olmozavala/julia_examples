# Examples at: https://github.com/FluxML/model-zoo/
using Flux

## Basics gradients
println("======== Basics: Gradients ========")
print("Gradient is 6x + 2 (14 for x = 2) Answer: ", gradient(x -> 3x^2 + 2x + 1, 2)) # 5 is the value of x
## Basics gradients
println("======== Basics: Gradients 2 ========")
f(x) = 3*x^2 + 2x + 1
df(x) = gradient(f, x)  # df/dx = 6x + 2

# println("Gradient of f(2) $(df(2))")
println("Gradient of f(2) $df")

## Multiple parameters
println("======== Basics: Multiple parameters ========")
f(x,y) = (3x - y)^2 # It returns a tuple with the gradient of each parameter (∂f/∂x₁, ∂f/∂x₂, ∂f/∂x₃,  .... ∂f/∂xₙ)
df(x,y) = gradient(f, x, y)

println("∇f (2(3x - y)*3 , 2(3x - y)*-1) ∇f(1,2) = ", df(1,2))

## Without passing the parameters
println("======== Basics: Multiple parameters ========")
x = [1,2] # [2 ,1] 

gs = gradient(Params(x,y)) do 
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