using DiffEqOperators, DifferentialEquations, OrdinaryDiffEq
using Plots

function getAdvectionModel()
    ## ====== Solving wave equation uₜₜ = c²uₓₓ
    Δx = .5
    x = Δx:Δx:2*π
    n = length(x)
    u₀ = sin.(x)
    du₀ = ones(n)

    order = 2 # approximation order
    deriv = 2 # Order of derivative

    ## ------- Define Boundary conditions -------
    bc = Dirichlet0BC(Float64)

    Δ = CenteredDifference(deriv, order, Δx, n)

    # Solve an ODE uₜₜ = c²uₓₓ
    # z = uₜ
    c = 2
    function wave(du,u,p,t) 
        # Because we need to add the variable z = du/dt, we need to add
        # another verctor with the same size of the spatial dimension of U 
        # In this case u' is on the frist n elements, and u on the other n
        du[1:n] = u[1:n]  # du[2] = zₜ = uₜₜ 
        du[n:n*2] = Δ*bc*(c^2 .* u[n:end])  # du[2] = zₜ = uₜₜ 
    end

    v₀ = [u₀... u₀...]
    tspan = (0.0, 10.0)
    prob = ODEProblem(wave, v₀, tspan)
    return init(prob, saveat=.1)
end