using Zygote
using Plots

## Gradients https://fluxml.ai/Zygote.jl/latest/#Taking-Gradients-1
# Basic use of gradient(f, x)
f(x,y) = x^2 + y^2
println("The gradient of f (x² + y²) =  is $(gradient(f,1, 1)))")

# gradient of dense network layer
linear(θ, x) = θ[:W]*x + θ[:b]
θ = Dict(:W => rand(4,2), :b => rand(4))
x = rand(2)
gradient(θ -> sum(linear(θ, x)), θ)[1]

## 
