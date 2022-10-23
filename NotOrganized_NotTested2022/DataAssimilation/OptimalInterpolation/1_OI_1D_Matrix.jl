
using LinearAlgebra
using Plots
using Random, Distributions
gr()

# This is an example of the Kalman Filter method for a 1D case. Here we assume
# that we can compute the observation ErrorCM and the background ECM
# Base.run(`clear`)
function plot_all(truth, obs, model, analysis)
    l = @layout [a b{0.20w} c d]
    p1 = heatmap(truth, c=:haline, title="Truth")
    p2 = heatmap(obs,   c=:haline, title="Obs")
    p3 = heatmap(model, c=:haline, title="Model")
    p4 = heatmap(model, c=:haline, title="Analysis")
    display(plot(p1,p2,p3,p4, layout=l, size=(1000,400)))
    xlabel!("Domain")
    ylabel!("Time")

end

function nantozero(x)
    if isnan(x)
        return 0
    else
        return x
    end
end

##  This code contains a 1D example of the Kalman Filter, following the notes from Kalnay  (page 179/369)
# Here we simulate that we only have observations on the first
# half of the domain. The other domain we don't

# 1D Kalman filter Tₐ = Tᵦ + W(T₀ - Tᵦ) 

# Time values
x = 0:.1:π  # This is the number of variables we have 
t = 0:.1:2π
tsteps = length(t)  # Simulated timesteps
n = length(x)       # Number of variables (model dimension)
p = Int(ceil(n/2))               # Number of observations (observations dimension)

# Simulated real state
Tₜ = zeros(tsteps, length(x))
for i=1:tsteps
    Tₜ[i,:] = sin.(x .+ t[i])
end
heatmap(x, t, Tₜ, label="True")
xaxis!("space")
yaxis!("time")

## Observations  (D1/2CD1/2
σₒ = ones(p)*.1  # Sigma for each variable
Cₒ = I(p) # Correlation matrix of erros (Independent errors)
Dₒ = Diagonal(σₒ.^2) # Covariance values at the diagonal
R = sqrt.(Dₒ)*Cₒ*sqrt.(Dₒ)# Rₚₓₚ Observation error covariance matrix of 
Tₒ = rand(MvNormal(σₒ), tsteps)' .+ Tₜ[:,1:p]  # Observation vector

# Background (Initial guess)
σᵦ = ones(n)*.1
Cᵦ = I(n) # Correlation matrix of erros (Independent errors)
Dᵦ = Diagonal(σᵦ.^2)
B = sqrt.(Dᵦ)*Cᵦ*sqrt.(Dᵦ) # Bₙₓₙ Background error covariance matrix of (model)
Tᵦ = rand(MvNormal(σᵦ), tsteps)' .+ Tₜ # Model vector

# Observation operator Hₚₓₙ (in this case it simply removes all the additional fiedls)
H = zeros(p,n)
H[1:p,1:p] = I(p)
Hᵀ = zeros(n,p) # Check this

# Optimal weight (nxp) constant
W = (B*H')*inv(H*B*H'+ R)
W = map(nantozero, W)
println("Optimal weight for σₒ=$(max(σₒ...)), σᵦ=$(max(σᵦ...)) is W=$(max(W...))")

# Analysis cycle (all together)
t=1
Tₐ = zeros(tsteps, n)
for t in 1:tsteps
    # Obtain current value of background (model) and observation
    c_Tᵦ = Tᵦ[t,:]
    c_Tₒ = Tₒ[t,:]
    Tₐ[t,:] = c_Tᵦ + W*(c_Tₒ - H*c_Tᵦ)
end
# Analysis error covariance
Pₐ = (I(n) - W*H)*B

## Plotting the results
plot_all(Tₜ, Tₒ, Tᵦ, Tₐ)

println("MSE Model: ", sum((Tₜ - Tᵦ).^2)/tsteps)
println("MSE Analysis: ", sum((Tₜ - Tₐ).^2)/tsteps)
println("Analysis variance: ", max(Pₐ...))