using ModelingToolkit, DiffEqOperators
using Plots
# Method of lines: http://www.scholarpedia.org/article/Method_of_lines

## Define Parameters, variables, and derivatives (ModelingToolkit)
@parameters t x  # Variables we know  (independent)
@variables u(..)
Dt = Differential(t)  # Derivative with respect to t
Dxx = Differential(x)^2  # Second derivative with respect to x

# 1D PDE and boundary conditions Δϕ = f -->  ∇²ϕ = f
eq  = Dt(u(t,x)) ~ Dxx(u(t,x))
bcs = [u(0,x) ~ cos(x),   
       u(t,0) ~ exp(-t),
       u(t,1) ~ exp(-t) * cos(1)]

# Space and time domains
domains = [t ∈ IntervalDomain(0.0,1.0),
        x ∈ IntervalDomain(0.0,1.0)]

# PDE system (https://mtk.sciml.ai/stable/systems/PDESystem/#)
pdesys = PDESystem(eq,bcs,domains,[t,x],[u])

# Method of lines discretization
dx = 0.1
order = 2
discretization = MOLFiniteDifference(dx,order)

# Convert the PDE problem into an ODE problem
prob = discretize(pdesys,discretization)

## Solve ODE problem
using OrdinaryDiffEq
sol = solve(prob,Tsit5(),saveat=0.2)

# Plot results and compare with exact solution
x = prob.space[2]
t = sol.t

# Method of Manufactured Solutions: exact solution
u_exact = (x,t) -> exp.(-t) * cos.(x)

plt = plot()
for i in 1:length(t)
    plot!(x,Array(prob.extrapolation[1](t[i])*sol.u[i]),label="Numerical, t=$(t[i])")
    scatter!(x, u_exact(x, t[i]),label="Exact, t=$(t[i])")
end
display(plt)
