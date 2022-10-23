
using LinearAlgebra
using Plots
using Random, Distributions
using Printf
using DiffEqBase

# This is an example of the Kalman Filter method for a 1D case for the Advection equation
# 1D Kalman filter Tₐ = Tᵦ + W(T₀ - Tᵦ) 

gr()
function plotCurrent(x, Tₐ, Tₜ, Tᵦ, Tₒ, t, Mfree)
    println("In plot: x=$(size(x)) Tₐ=$(size(Tₐ)) Tₜ=$(size(Tₜ)) Tᵦ=$(size(Tᵦ)) Tₒ=$(size(Tₒ)) Mfree=$(size(Mfree))")
    # Simple helper function to plot all the states together
    plot(x, Tₜ, label="True", c=:gray, title="Time $(round(t,digits=3))", ylim=(-1.1,1.1), size=(1000,1000))
    scatter!(x, Tᵦ, label="Background", c=:blue, ms=3, markerstrokewidth=0)
    scatter!(x, Tₒ, label="Observations", c=:cyan, ms=3, markerstrokewidth=0)
    plot!(x, Mfree, label="FreeModel", c=:red)
    display(plot!(x, Tₐ, label="Analysis", c=:green))
end
## 

# Time values
x = 0:.1:2*π  # This is the number of variables we have 
tend = 100.0
Δt = 0.1
tsteps = Int(tend/Δt)
n = length(x)     # Number of variables (model dimension)
p = Int(n/1)      # Number of observations (observations dimension)
N = 50            # Number of ensembles

# --- Making the 'true state' ---
v = 2 # Velocity of the advection
Tₜall = zeros(tsteps, length(x))
for i=1:tsteps
    Tₜall[i,:] = sin.(x .+ v*Δt*(i-1))
end
# heatmap(x, range(0.0, tend, length=tsteps), Tₜ, title="True field",
                # label="True", xaxis="x", yaxis="time")
plot(x, Tₜall[1,:], title="True state at t=0")

## ================ Running KF =====================

# ---- Observations
H = I(p)# Observation operator Hₚₓₙ
σₒ = ones(p)*.05  # Standard deviation for each variable (same for all)
Cₒ = I(p) # Correlation matrix of errors (Independent errors)
Dₒ = Diagonal(σₒ.^2) # Variance of the error 
R = sqrt.(Dₒ)*Cₒ*sqrt.(Dₒ)# Rₚₓₚ Observation error covariance matrix (static through time)

##
include("../AdvectionModel.jl")
Mfree = getAdvectionModel(x, 0.0, tend, Δt, v) # Get model
σperturb = .1
Mensemble = [getAdvectionModel(x, 0.0, tend, Δt, v, σperturb) for i in 1:N]

for tᵢ=1:tsteps
# for tᵢ=1:20
    t =Mensemble[1].t

    # ======================= ANALYSIS STEP ==============
    Tₜ = Tₜall[tᵢ,:]
    Tₒ = rand(MvNormal(σₒ), 1) .+ H*Tₜ  # Observation vector

    # ---------- Computing Background cov err from the MSE of the ensemble ----
    Tᵦ = reduce(+, [M.u for M in Mensemble])/N  # Mean state of the ensemble
    B = (1/(N+1))*reduce(+,[(M.u - Tᵦ)*(M.u - Tᵦ)' for M in Mensemble])
    if tᵢ == 1
        display(heatmap(B))
    end

    # ----------- Make the analysis -----
    K = (B*H')*inv(H*B*H'+ R)
    Tₐ = Tᵦ + K*(Tₒ - H*Tᵦ)
    Pₐ = (I(n) - K*H)*B

    # -----------   Plot results ----------
    if (tᵢ % 10 == 0) || tᵢ == 1
        println(" -------- Analysis at time = $(round(t,digits=3)) ---------")
        plotCurrent(x, Tₐ, Tₜ, Tᵦ, H*Tₒ, t, Mfree.u)
        println("MSE Free Model: ", round(sum((Tₜ - Mfree.u).^2)/n, digits=4))
        println("MSE Analysis EKF: ", round(sum((Tₜ - Tₐ).^2)/n, digits=4))
        # println("Analysis variance: ", max(Pₐ...))
        sleep(.3)
    end
    # Replace current U with the analysis one

    # -----------   Update each ensemble ----------
    for i=1:N
        DiffEqBase.set_u!(Mensemble[i],Tₐ[:,1] .+ (rand(n).-.5).*σperturb)
        step!(Mensemble[i], Δt, true)
    end
    step!(Mfree, Δt, true)
end