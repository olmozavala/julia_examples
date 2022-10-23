using DiffEqOperators, DifferentialEquations, OrdinaryDiffEq
using Plots

## ========= Advection 1D  ∂f/∂t + ∇ ⋅ (fu) = 0  → ∂f/∂t + (uᵢfₓ + uⱼf𝒚 + uₖf𝒛) = 0  
function getAdvectionModel(x, t0, t1, Δt, v, Δσ=0)
    ### x -> x Axis, t0 and t1 time range, v velocity, Δσ if we want to perturb the initial state

    n = length(x)
    Δx = x[2]-x[1]
    u₀ = sin.(x) + (rand(n).-.5).*Δσ # This is our scalar field
    ##
    ord_deriv = 1
    ord_approx = 2
    # Δuw = UpwindDifference(ord_deriv, ord_approx, Δx, n, 1)
    Δuw = CenteredDifference(ord_deriv, ord_approx, Δx, n)

    U = v*ones(n) # Vector field moving homogeneously to one side
    bcx = PeriodicBC(Float64)

    stepuw(u,p,t) = Δuw*bcx*(u.*U)  # Define the ODE with UpwindDifference

    # Solving equation
    probuw = ODEProblem(stepuw, u₀, (t0, t1))
    # solve(probuw, saveat=Δt)
    return init(probuw, Tsit5(), saveat=Δt)
end

## Only for testing
# n = 10
# Δx = 2*π/n
# Δt = .5
# x = range(Δx, step=Δx, length=n) 

# M = getAdvectionModel(x, 0.0, 20.0, Δt, 2) # Get model

# for i=1:100
#     println("Time= $(M.t) U=$(M.u[1])")
#     display(plot(x, M.u, title="T = $(M.t)"))
#     step!(M, Δt, true)
# end