using LinearAlgebra
using Plots
using Random, Distributions
gr()
Base.run(`clear`)

## 

# Single particle Kalman filter Tₐ = Tᵦ + W(T₀ - T\_b) 

# Time values
x = 0:.1:2π
tsteps = length(x)

# Simulated real state
Tₜ = sin.(x)  # Simulated real state

# Observations
σₒ = .2
Tₒ = rand(Normal(0,σₒ), tsteps) .+ Tₜ

# Background (Initial guess)
σᵦ = .2
Tᵦ = rand(Normal(0,σᵦ), tsteps) .+ Tₜ

# Optimal weight
W = σᵦ.^2/(σᵦ.^2 + σₒ.^2)
println("Optimal weight for σₒ=$σₒ, σᵦ=$σᵦ is W=$W")

# Analysis cycle (all together)
Tₐ = Tᵦ + W*(Tₒ - Tᵦ)

##
p = plot(x, Tₜ, label="True", color = :cyan, xaxis="Time", yaxis="u")
scatter!(x, Tₒ, label="Obs", msw = 0, mc = :blue, ms=2)
scatter!(x, Tₐ, label="Analysis", msw = 0, mc = :green, ms=4)
scatter!(x, Tᵦ, label="Model", msw = 0, mc = :red, ms=2)
display(p)

println("MSE Model: ", sum((Tₜ - Tᵦ).^2)/tsteps)
println("MSE Analysis: ", sum((Tₜ - Tₐ).^2)/tsteps)