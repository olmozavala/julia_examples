## 
using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DomainSets, Plots, Printf
gr()

## --------- MethodOfLines Use
# https://docs.juliahub.com/MethodOfLines/v6ocr/0.3.0/MOLFiniteDifference/#molfd

@parameters x t
dx = 0.1
println("Discretizing problem...")
# MOLFiniteDifference(dxs, time, approx_order::Int, upwind_order::Int, grid_align::G)
discretization = MOLFiniteDifference([x=>dx],t, approx_order = 2)
println(discretization)
##

## --------- Heat equation Example --------
@parameters t, x
@variables u(..)

Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2  # Second derivative

dx = 0.1
α = 1.1

eq = [Dt(u(t,x)) ~ α * Dxx(u(t,x))]

domain = [x ∈ Interval(0.0, 10.0),
          t ∈ Interval(0.0, 1.0)]

# --- Plotting intial conditions
ic(x) = exp(-(x - 4.0)^2) +  2*exp(-(x - 8.0)^2)
a = 0:0.1:10
p = plot(a, ic.(a), title="Initial conditions", label="u(x,0)", legend=:top)
display(p)

ic_bc = [u(0, x) ~ exp(-(x - 4.0)^2) +  2*exp(-(x - 8.0)^2),
         u(t, 0.0) ~ 0.0,
         u(t, 10.0) ~ 0.0]

@named sys = PDESystem(eq,ic_bc,domain,[t,x],[u(t,x)])

##
print("Discretizing problem...")
dx = 0.1
discretization = MOLFiniteDifference([x=>dx],t, approx_order = 2)
prob = discretize(sys, discretization)
print("Done!")
##
print("Solving ODE...")
Δt = 0.1
sol = solve(prob, Tsit5(), saveat=0.01) # Specify the error tolerance (compare and see)
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
    xaxis="X", yaxis="Time", ylim=(0,4), 
    title=title, label="u(x,t)")
end
## As heatmap
heatmap(xdom, tdom, sol.u[u(t,x)], xaxis="X", yaxis="Time")