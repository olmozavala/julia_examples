# Examples at: https://fluxml.ai/Zygote.jl/latest/
using Zygote
using Plots

## Basics gradients inline
println("======== Basics: Gradients inline ========")
print("Gradient is 6x + 2 (14 for x = 2) Answer: ", gradient(x -> 3x^2 + 2x + 1, 2)) # 5 is the value of x

## Basics gradients from function
println("======== Basics: Gradients by variable ========")
f(x) = 3*x^2 + 2x + 1
df(x) = gradient(f, x)  # df/dx = 6x + 2
println("Gradient df = 6x + 2 f(2) $(df(2))")

## Multiple parameters
println("======== Basics: Multiple parameters ========")
f(x,y) = (3x - y)^2 # It returns a tuple with the gradient of each parameter (∂f/∂x₁, ∂f/∂x₂, ∂f/∂x₃,  .... ∂f/∂xₙ)
df(x,y) = gradient(f, x, y)
println("∇f = 6(3x - y) , -2(3x - y) ∇f(1,2) = ", df(1,2))

##
println("======== Gradients for arrays (as long the function returns a scalar)========")
f(x) = sum(3x.^2)
df(x) = gradient(f, x)
println("∇f = (6x₁, 6x₂, 6x₃, ...) = ∇f([3,5,8])= ", df([3,5,8]))

## Linear network 
println("======== Gradients example for linear networks========")
linear(θ, x) = θ[:W]*x .+ θ[:b]

in_vars = 5
x = rand(in_vars)
θ = Dict(:W => rand(1,in_vars), :b => rand(in_vars))
# Gradient of θ with respect to x
gradient(θ -> sum(linear(θ, x)), θ)

## Implicit parameters Flux style
println("======== Implicit parameters Flux type =========")
x = rand(3)
W = rand(1, 3); b = rand(1);
linear(x) = W * x .+ b
grads = gradient(() -> sum(linear(x)), Params([W, b]))
grads[W], grads[b]