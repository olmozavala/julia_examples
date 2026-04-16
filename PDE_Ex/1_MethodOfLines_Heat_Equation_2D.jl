## --------- 2D Heat equation Example --------
using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DomainSets, Plots, Printf
gr()

## --------- MethodOfLines Use
# https://docs.juliahub.com/MethodOfLines/v6ocr/0.3.0/MOLFiniteDifference/#molfd

# 1. Define the symbolic parameters x, y, and t for the PDE problem.
# 2. Set the spatial discretization step sizes dx and dy.
# 3. Set up the discretization scheme using MethodOfLines: 
#    - spatial steps [x=>dx, y=>dy], time variable t, 2nd order finite difference.

## --------- PDE System ---------
@parameters t, x, y
@variables u(..)

Dt = Differential(t)
Dx = Differential(x)
Dy = Differential(y)
Dxx = Differential(x)^2  # Second derivative in x
Dyy = Differential(y)^2  # Second derivative in y

dx = 0.2
dy = 0.2
α = 1.1 # Diffusivity coefficient

# Define the 2D heat equation: uₜ = α * (uₓₓ + u_yy)
eq = [Dt(u(t,x,y)) ~ α * (Dxx(u(t,x,y)) + Dyy(u(t,x,y)))]

domain = [x ∈ Interval(0.0, 10.0),
          y ∈ Interval(0.0, 10.0),
          t ∈ Interval(0.0, 1.0)]

## --- Initial and Boundary conditions
# Initial condition: Two 2D Gaussian peaks
ic_func(x, y) = exp(-((x - 3.0)^2 + (y - 3.0)^2)) + 1.5*exp(-((x - 7.0)^2 + (y - 7.0)^2))

ic_bc = [u(0, x, y) ~ ic_func(x, y),
         u(t, 0.0, y) ~ 0.0,      # Left boundary
         u(t, 10.0, y) ~ 0.0,     # Right boundary
         u(t, x, 0.0) ~ 0.0,      # Bottom boundary
         u(t, x, 10.0) ~ 0.0]     # Top boundary

## --- Plotting initial conditions
x_grid = 0:0.2:10
y_grid = 0:0.2:10
ic_values = [ic_func(xi, yi) for yi in y_grid, xi in x_grid]

p = plot(layout = (1,2), size=(900,400))
# Initial condition as heatmap
heatmap!(p[1], x_grid, y_grid, ic_values, 
         title="Initial Condition: u(x,y,0)", 
         xlabel="x", ylabel="y", 
         c=:viridis)
# Initial condition as surface
surface!(p[2], x_grid, y_grid, ic_values, 
         title="Initial Condition: u(x,y,0)", 
         xlabel="x", ylabel="y", zlabel="u",
         c=:viridis)
display(p)

## --------- PDE System and discretization---------
@named sys = PDESystem(eq,ic_bc,domain,[t,x,y],[u(t,x,y)])

print("Discretizing problem...")
discretization = MOLFiniteDifference([x=>dx, y=>dy], t, approx_order = 2)
prob = discretize(sys, discretization)
println("Done!")

## --------- Solving the PDE problem---------
print("Solving ODE...")
sol = solve(prob, Tsit5(), saveat=0.02) # Specify the save time interval
println("Done!")

## ----- Analyzing the solution object  https://github.com/SciML/SciMLBase.jl/blob/master/src/solutions/pde_solutions.jl
println("List of independent variables: ", sol.ivs)
println("List of dependent variables: ", sol.dvs)
println("Domains: ", sol.ivdomain)
println("Number of time values: ",  length(sol.t))
tdom = sol.ivdomain[1]
xdom = sol.ivdomain[2]
ydom = sol.ivdomain[3]

## First time step - heatmap
heatmap(xdom, ydom, sol.u[u(t,x,y)][1,:,:], 
        xaxis="X", yaxis="Y", 
        title="u(x,y,0)", c=:viridis)

## First time step - surface plot
surface(xdom, ydom, sol.u[u(t,x,y)][1,:,:], 
        xaxis="X", yaxis="Y", zlabel="u",
        title="u(x,y,0)", c=:viridis)

## Animation as heatmap
@gif for i ∈ 1:length(tdom)
    title_text = @sprintf "Time %0.2f" tdom[i]
    heatmap(xdom, ydom, sol.u[u(t,x,y)][i,:,:], 
            xaxis="X", yaxis="Y", 
            title=title_text, 
            clim=(0, 2.0),
            c=:viridis)
end

## Animation as surface plot
@gif for i ∈ 1:length(tdom)
    title_text = @sprintf "Time %0.2f" tdom[i]
    surface(xdom, ydom, sol.u[u(t,x,y)][i,:,:], 
            xaxis="X", yaxis="Y", zlabel="u",
            title=title_text, 
            zlim=(0, 2.0),
            c=:viridis)
end