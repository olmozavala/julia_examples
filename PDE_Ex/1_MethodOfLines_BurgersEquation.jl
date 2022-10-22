## 
using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DomainSets, Plots, Printf
gr()

## --------- 
@parameters t, x
@variables u(..)

Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2  # Second derivative


eq = [Dt(u(t,x)) ~ -u(t,x)*Dx(u(t,x)) + (0.01/π)*Dxx(u(t,x))]

domain = [x ∈ Interval(-1.0, 1.0),
          t ∈ Interval(0.0, 1.0)]


ic_bc = [u(0, x) ~ -sin(π*x),
         u(t, -1.0) ~ 0.0,
         u(t, 1.0) ~ 0.0]

@named sys = PDESystem(eq,ic_bc,domain,[t,x],[u(t,x)])

##
print("Discretizing problem...")
dx = 0.005
discretization = MOLFiniteDifference([x=>dx],t, approx_order = 2)
prob = discretize(sys, discretization)
print("Done!")

##
print("Solving ODE...")
Δt = 0.05
sol = solve(prob, Tsit5(), saveat = 0.01) # Specify the error tolerance (compare and see)
print("Done!")

## ----- Analyzing the solution object  https://github.com/SciML/SciMLBase.jl/blob/master/src/solutions/pde_solutions.jl
println("List of independent variables:", sol.ivs)
println("List of dependent variables:", sol.dvs)
println("Domains:", sol.ivdomain)
println("Time values",  sol.t)
tdom = sol.ivdomain[1]
xdom = sol.ivdomain[2]

## First time step
plot(xdom, sol.u[u(t,x)][1,:], xaxis="X", yaxis="Time")
## All time steps
# As animation
@gif for i ∈ 1:length(tdom)
    title = @sprintf "Time %0.2f" tdom[i]
    p = plot(xdom, sol.u[u(t,x)][i,:], 
    xaxis="X", yaxis="Time",
    title=title, label="u(x,t)", ylim=[-1,1])
end
## As heatmap
# heatmap(xdom, tdom, sol.u[u(t,x)], xaxis="X", yaxis="Time")
heatmap(tdom, xdom, transpose(sol.u[u(t,x)]), xaxis="Time", yaxis="X")