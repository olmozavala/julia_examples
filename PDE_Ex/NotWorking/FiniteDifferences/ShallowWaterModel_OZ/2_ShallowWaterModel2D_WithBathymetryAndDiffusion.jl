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
h_topography = repeat(range(0.0,0.5,length=n)', n, 1) 

# Boundary conditions  au + b ∂u/∂n = c  
dirx = 0.0
diry = 0.0
dir_bcx = RobinBC{1}((1.0, 0.0, dirx), (1.0, 0.0, diry), (Δx), 1, (n,n)) # Defines Dirichlet BC with 0
dir_bcy = RobinBC{2}((1.0, 0.0, dirx), (1.0, 0.0, diry), (Δy), 1, (n,n)) # Defines Dirichlet BC with 0
newman_bcx = RobinBC{1}((0.0, 1.0, 0.0), (0.0, 1.0, 0.0), (Δx), 1, (n,n)) # Defines Neumann = 0 BC x 
newman_bcy = RobinBC{2}((0.0, 1.0, 0.0), (0.0, 1.0, 0.0), (Δy), 1, (n,n)) # Defines Neumann = 0 BC y

# heatmap(x, y, h₀, label="SSH",  markersize=7, markerstrokewidth=0)
a = plot(x, y, h₀, title="SSH",  c=:blues, st=:surface, camera=(20,60), alpha=0.5)
b = plot(x, y, h_topography, title="Topography",  c=:greens, st=:surface, camera=(20,60), alpha=0.5)
c = scatter([x[1]-Δx x... x[end]+Δx]', (newman_bcx*h₀)[:,Int(round(n/2))], label="SSH BC", title="SSH BC",  c=:blue)
plot(a,b, layout=@layout[a b], size=(900, 300))

## We set f in the following way [u h]
ord_deriv = 1
ord_approx = 2
Δxop = CenteredDifference{1}(1, ord_approx, Δx, n)
Δyop = CenteredDifference{2}(1, ord_approx, Δy, n)
Δxdiff = CenteredDifference{1}(2, 2, Δx, n)
Δydiff = CenteredDifference{2}(2, 2, Δy, n)

# We organize f in the following way u, v, h
f₀ = [zeros(n,n); zeros(n,n); h₀];
##
function step(du,f,p,t)
    n = p[1]
    u = @view f[1:n, :]
    v = @view f[n+1:2*n, :]
    ssh = @view f[2*n+1:3*n, :]
    g = 9.81
    n = p[1]
    H = 1
    du[1:n      , :] = -g*Δxop*newman_bcx*(ssh .- h_topography) + Δxdiff*dir_bcx*u # Solving for u
    du[n+1:2*n  , :] = -g*Δyop*newman_bcy*(ssh .- h_topography) + Δydiff*dir_bcy*v   # Solving for v
    # du[1:n      , :] = -g*Δxop*newman_bcx*(ssh .- h_topography) # Solving for u
    # du[n+1:2*n  , :] = -g*Δyop*newman_bcy*(ssh .- h_topography) # Solving for v
    du[2*n+1:3*n, :] = -H*(Δxop*dir_bcx*u +  Δyop*dir_bcy*v) # Solving for h
end

# Solving equation
t0 = 0.0
t1 = 0.2
Δt = .01
tsteps = Int(round(t1/Δt)) - 1
println("Solving with t=($t0, $t1) Δt=$Δt and Δx=$Δx ....")
p = (n) # Parameters (just the size of n)
prob = ODEProblem(step, f₀, (t0, t1), p)
sol = solve(prob, Tsit5(), saveat=Δt)
println("Done!")

##  
println("Plotting...")
lw = 4
mid_i =  Int(round(n/2))#Middle index
for i=1:Int(round(tsteps/10)):tsteps
# for i=1:10
    a = heatmap(x, y, sol.u[i][1:n, :], title="U t = $(sol.t[i])")
    b = heatmap(x, y, sol.u[i][n+1:2*n, :],  title="V t = $(sol.t[i])")
    c = plot(x, y, sol.u[i][2*n+1:3*n, :], title="SSH t = $(sol.t[i])",  
                        c=:blues, st=:surface, camera=(20,60), alpha=0.5)
    d = plot(x, sol.u[i][2*n+1:3*n, mid_i] .+ h_topography[mid_i,:], title="SSH y=$mid_i t = $(sol.t[i])")
    p = plot(a,b,c,d, layout=@layout[a b; c d], size=(700,700))
    display(p)
    sleep(.01)
end

println("Done")