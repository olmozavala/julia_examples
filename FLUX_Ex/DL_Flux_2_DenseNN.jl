# Examples at: https://github.com/FluxML/model-zoo/
using Flux
using Plots
gr() # Using GR to plot

## ------------- Building 'Dense' layer manually
# https://fluxml.ai/Flux.jl/stable/models/basics/#Stacking-It-Up-1
struct Affine
  W
  b
end

Affine(in::Integer, out::Integer) =
  Affine(randn(out, in), randn(out))

# Overload call, so the object can be used as a function
(m::Affine)(x) = m.W * x .+ m.b

a = Affine(10, 5)

a(rand(10)) # => 5-element vector
print(a)

## ------------- Using Flux Dense layers

# Side note, sigmoid function can be called with σ
plot(σ.(-10:10))

# Chain Example
m = Chain(x-> x^2, x -> x+1)
print(m(3))

# x = rand(10)
# scatter(x, label="original")
# scatter!(softmax(x), label="softmax")

model1 = Chain(
    Dense(10,5, σ),
    Dense(5,1),
    softmax
    )
    
#softmax is also a function
# plot(softmax(-10:10))

println(model1(rand(10)))