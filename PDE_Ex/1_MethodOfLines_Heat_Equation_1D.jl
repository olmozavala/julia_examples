## --------- 1D Heat equation Example --------
using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DomainSets, Plots, Printf
gr()

## --------- MethodOfLines Use
# https://docs.juliahub.com/MethodOfLines/v6ocr/0.3.0/MOLFiniteDifference/#molfd

# 1. Define the symbolic parameters x and t for the PDE problem.
# 2. Set the spatial discretization step size dx to 0.1.
# 3. Set up the discretization scheme using MethodOfLines: 
#    - spatial step [x=>dx], time variable t, 2nd order finite difference.

# --------- Discretization ---------
@parameters x t
dx = 0.1
println("Discretizing problem...")
# MOLFiniteDifference(dxs, time, approx_order::Int, upwind_order::Int, grid_align::G)
discretization = MOLFiniteDifference([x=>dx],t, approx_order = 2)
println("Discretization: ", discretization)

## --------- PDE System ---------
@parameters t, x
@variables u(..)

Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2  # Second derivative

dx = 0.1
α = 1.1 # Diffusivity coefficient

# 4. Define the PDE equation: uₜ = α * uₓₓ
eq = [Dt(u(t,x)) ~ α * Dxx(u(t,x))]

domain = [x ∈ Interval(0.0, 10.0),
          t ∈ Interval(0.0, 1.0)]

## --- Plotting intial conditions
ic(x) = exp(-(x - 4.0)^2) +  2*exp(-(x - 8.0)^2)
a = 0:0.1:10

ic_bc = [u(0, x) ~ exp(-(x - 4.0)^2) +  2*exp(-(x - 8.0)^2),
         u(t, 0.0) ~ 0.0,
         u(t, 10.0) ~ 0.0]

p = plot(layout = (3,1), size=(400,900))  # change to a column layout

# Initial condition at t=0 (center plot)
plot!(p[1], a, ic.(a), title="Initial Condition: u(x,0)", xlabel="x", ylabel="u", label="u(x,0)", legend=:top)
# Boundary condition at x=0: u(t,0) = 0 (left plot)
plot!(p[2], [0, 0], [0, 1], title="Boundary: u(t,0)", xlabel="x", ylabel="t", label="u(t,0)=0", legend=:top)
# Boundary condition at x=10: u(t,10) = 0 (right plot)
plot!(p[3], [10, 10], [0, 1.0], xlim=(0,10), ylim=(0,1), title="Boundary: u(t,10)", xlabel="x", ylabel="t", label="u(t,10)=0", legend=:top)

display(p)
## --------- PDE System and discretization---------
@named sys = PDESystem(eq,ic_bc,domain,[t,x],[u(t,x)])

print("Discretizing problem...")
dx = 0.1
discretization = MOLFiniteDifference([x=>dx],t, approx_order = 2)
prob = discretize(sys, discretization)
print("Done!")

## --------- Solving the PDE problem---------
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
# As animation
@gif for i ∈ 1:length(tdom)
    title = @sprintf "Time %0.2f" tdom[i]
    p = plot(xdom, sol.u[u(t,x)][i,:], 
    xaxis="X", yaxis="Time", ylim=(0,4), 
    title=title, label="u(x,t)")
end
## As heatmap
heatmap(xdom, tdom, sol.u[u(t,x)], xaxis="X", yaxis="Time")