using ModelingToolkit, DiffEqOperators, DifferentialEquations
using Plots
using Random
using DiffEqOperators, OrdinaryDiffEq
gr()
include("../../../OZ_Tools/ExampleData.jl")

## ========= Simple ShallowWaterModel 2D ================
n = 50# Number of spatial grid points
Δx = 2*π/n
Δy = 2*π/n
x = range(Δx, step=Δx, length=n) 
y = range(Δy, step=Δy, length=n) 
h₀ = getGaussShape(n, n/20, 2)

# Boundary conditions  au + b ∂u/∂n = c  
dirx = 0.0
diry = 0.0
dir_bcx = RobinBC{1}((1.0, 0.0, dirx), (1.0, 0.0, diry), (Δx), 1, (n,n)) # Defines Dirichlet BC with 0
dir_bcy = RobinBC{2}((1.0, 0.0, dirx), (1.0, 0.0, diry), (Δy), 1, (n,n)) # Defines Dirichlet BC with 0`
newman_bcx = RobinBC{1}((0.0, 1.0, 0.0), (0.0, 1.0, 0.0), (Δx), 1, (n,n)) # Defines Neumann = 0 BC x 
newman_bcy = RobinBC{2}((0.0, 1.0, 0.0), (0.0, 1.0, 0.0), (Δy), 1, (n,n)) # Defines Neumann = 0 BC y

# heatmap(x, y, h₀, label="SSH",  markersize=7, markerstrokewidth=0)
a = plot(x, y, h₀, title="SSH",  c=:blues, st=:surface, camera=(20,60), alpha=0.5)
b = scatter([x[1]-Δx x... x[end]+Δx]', (newman_bcx*h₀)[:,Int(round(n/2))], label="SSH BC", title="SSH BC",  c=:blue)
plot(a,b, layout=@layout[a b], size=(800, 400))

## We set f in the following way [u h]
ord_deriv = 1
ord_approx = 2
Δxop = CenteredDifference{1}(ord_deriv, ord_approx, Δx, n)
Δyop = CenteredDifference{2}(ord_deriv, ord_approx, Δy, n)

# We organize f in the following way u, v, h
f₀ = [zeros(n,n); zeros(n,n); h₀];
##
function step(du,f,p,t)
    n = p[1]
    u = f[1:n, :]
    v = f[n+1:2*n, :]
    ssh = f[2*n+1:3*n, :]
    g = 9.81 
    H = 1
    du[1:n      , :] = -g*Δxop*newman_bcx*ssh # Solving for u
    du[n+1:2*n  , :] = -g*Δyop*newman_bcy*ssh # Solving for v
    du[2*n+1:3*n, :] = -H*(Δxop*dir_bcx*u +  Δyop*dir_bcy*v) # Solving for h

    # With "manual BC" to simulate a Gulf of Mexico
    # Small open on the south east with 'input' and small output on the northeast with 'output'
    # n = 50
    # nmid = 25

    # du[2*n+1:3*n, :] = -H*(Δxop*[zeros(n) u' zeros(n)]' +  Δyop*dir_bcy*v) # Solving for h
end

# Solving equation
t0 = 0.0
t1 = 5.0
Δt = 0.001
tsteps = Int(round(t1/Δt)) - 1
println("Solving with Δt=$Δt and Δx=$Δx ....")
p = (n) # Parameters (just the size of n)
prob = ODEProblem(step, f₀, (t0, t1), p)
sol = solve(prob, saveat=Δt);
println("Done!")

##  
println("Plotting...")
lw = 4
for i=1:Int(round(tsteps/30)):tsteps
# for i=1:1023232151232
    # a = plot(x, y, sol.u[i][1:n, :], c=:blues, st=:surface, title="U t = $(sol.t[i])")
    # b = plot(x, y, sol.u[i][n+1:2*n, :], c=:blues, st=:surface, title="V t = $(sol.t[i])")
    # c = plot(x, y, sol.u[i][2*n+1:3*n, :], c=:blues, st=:surface, title="SSH t = $(sol.t[i])")
    a = heatmap(x, y, sol.u[i][1:n, :], title="U t = $(sol.t[i])")
    b = heatmap(x, y, sol.u[i][n+1:2*n, :],  title="V t = $(sol.t[i])")
    # c = heatmap(x, y, sol.u[i][2*n+1:3*n, :],  title="SSH t = $(sol.t[i])")
    c = plot(x, y, sol.u[i][2*n+1:3*n, :], title="SSH t = $(sol.t[i])",  c=:blues, st=:surface,
                        zlim=[-1.0,1.0], camera=(20,60), alpha=0.5)
    d = plot(x, sol.u[i][2*n+1:3*n, Int(round(n/2))], title="SSH t = $(sol.t[i])", ylim=[-1,1])
    p = plot(a,b,c,d, layout=@layout[a b; c d], size=(900,300))
    display(p)
    sleep(.01)
end

println("Done")



###o
dir_bcx = RobinBC{1}((1.0, 0.0, dirx), (1.0, 0.0, diry), (Δx), 1, (2,2)) # Defines Dirichlet BC with 0