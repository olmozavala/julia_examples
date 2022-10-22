# Examples at: https://github.com/FluxML/model-zoo/
using Flux
using Plots
gr() # Using GR to plot

## ------------- Building 'Dense' layer manually
# https://fluxml.ai/Flux.jl/stable/models/basics/#Stacking-It-Up-1
println("Dense layer manually")
struct Affine
  W
  b
end

# Defines an overload. When Affine is called with integers it will 
# create a struct with random matrix W and vector b
Affine(in::Integer, out::Integer) = Affine(randn(out, in), randn(out))

# Defines a function overload. It says, if you call an Affine structure with (x)
# then it will compute the linear function m.W * x .+ m.b
(m::Affine)(x) = m.W * x .+ m.b

a = Affine(10, 5) # Initialize your affine transformation randomly with 10 inputs and 5 outputs

println("The output of my affine transformation for 10 random numbers is:\n $(a(rand(10)))")

## ------------- Using Flux Dense layers
println("Built in dense layer")

# Side note, sigmoid function can be called with σ
plot(σ.(-10:10), titile="Sigmoid with σ")

# Chain Example
m = Chain(x-> x^2, x -> x+1)
println("Example use of chain, it chains multiple operations from left to right: $(m(3))")

# scatter(softmax(rand(10)), label="softmax")

model1 = Chain(
    Dense(10,5, σ), # Dense NN with 10 inputs, 5 outputs and sigmoid activation function
    Dense(5,1),  # Dense NN with 5 inputs one output linear activation function
    softmax  # Softmax transform the output into probability function
    )
    
println("Output of initial model with model1(rand(10)))")