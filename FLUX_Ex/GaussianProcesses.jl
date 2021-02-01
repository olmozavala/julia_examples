using Plots
using LinearAlgebra
using Statistics
gr() # Using GR to plot
## Imports 

# Here we are simulating that each point in our data comes from
# a gaussian distribution with mean at that value and 0 std

function kernel(a, b)
    # Kernel 1
    m = length(a)
    n = length(b)
    K = ones(m+n, m+n)
    for i in 1:m
        for j in 1:n
            K[i,j] = exp.( (a[i] - b[j])^2 )
        end
    end
    return K
end


## Super basics just 3 points
x = [0 1 2]
y = [0 2 4]
xₛ = 1.5

scatter(x,y)

Kall = kernel(x, xₛ)
Kss = kernel(xₛ, xₛ)
Ksᵀ = transpose(Ks)
print(Ks)

print("Size of K: $(size(K)))")
L = cholesky(K)

# μₛ = xₛ .+ Ksᵀ*K¹*(f .- μ)
# Σ = Kss - Ksᵀ*K¹*Ks

# scatter!(xₛ, μₛ, c="green")


## A little bit more real
examples = 10
x = -2π .+ rand(examples).*4π
f = sin.(x) 
μ = x
K = kernel(x, x)
scatter(x, f)

xₛ = -1π .+ rand(examples).*2π
scatter!(xₛ, 0 .* xₛ, color="red")
K¹ = inv(K)

Ks = kernel(x, xₛ)
Kss = kernel(xₛ, xₛ)
Ksᵀ = transpose(Ks)

print("Size of K: $(size(K)))")
# L = cholesky(K)

μₛ = xₛ .+ Ksᵀ*K¹*(f .- μ)
Σ = Kss - Ksᵀ*K¹*Ks

scatter!(xₛ, μₛ, c="green")