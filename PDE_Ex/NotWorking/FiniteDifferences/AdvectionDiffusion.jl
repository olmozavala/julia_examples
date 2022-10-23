
using ModelingToolkit, DiffEqOperators, DifferentialEquations
using Plots
using Random
using DiffEqOperators, OrdinaryDiffEq
gr()

## ========= Advection-Diffusion or Convection-Diffusion 1D  ∂f/∂t + ∇ ⋅ (fu) = 0  → ∂f/∂t + (uᵢfₓ + uⱼf𝒚 + uₖf𝒛) = 0  
n = 100# Number of spatial grid points
Δx = 2*π/n
x = range(Δx, step=Δx, length=n) 
u₀ = getGaussShape(n, n/8)
scatter(u₀)

##
order_approx = 2
Δdiff = CenteredDifference(2, ord_approx, Δx, n)
Δadv = UpwindDifference(1, ord_approx, Δx, n, 1)
display(Matrix(Δuw))
# println(size(Matrix(Δuw)))

v = 1
U = v*ones(n) # Vector field moving homogeneously to one side
# U = -cos.(x)# Vector field moving to both sides from the middle

# Boundary conditions  au + b ∂u/∂n = c  
# bcx = RobinBC((1.0, 0.0, 0.0),(1.0, 0.0, 0.0), (Δx), 1, size(u₀)) # Direichlet = 0
# bcx = Neumann0BC(Δx)
# bcx = Dirichlet0BC(Float64)
bcx = PeriodicBC(Float64)

scatter(x, u₀, label="Scalar field")
plot!(x, U, c=:green, label="Vector field", legend=:bottomright)
display(scatter!(transpose([x[1]-Δx x... x[end]+Δx]), bcx*u₀, label="BC", markersize=2, markerstrokewidth=0, mcolor=:red, c=:red))
step(u,p,t) = Δdiff*bcx*(u.*U) + Δadv*bcx*u

# Solving equation
t0 = 0.0
t1 = 2
Δt = .1
tsteps = Int(round(t1/Δt)) - 1
println("Solving with Δt=$Δt and Δx=$Δx")
prob = ODEProblem(step, u₀, (t0, t1))
sol = solve(prob, saveat=Δt);
println("Done!")

##
println("Plotting...")
plot(u₀, title="Solution t=0", ylim=[-.1, 1.1])
lw = 4
for i=1:Int(round(tsteps/10)):tsteps
    # p = plot(-sin.(x .+ v*(Δt*(i-1))), title="Solution True t=$(soluw.t[i])", label="True", c=:green, lw=lw)
    display(plot(sol[:,i], c=:blue, label="CD", ylim=[-.1, 1.1], title="Time $(sol.t[i])"))
    sleep(.1)
end
println("Done")


## ========= Advection 2D  ∂f/∂t + ∇ ⋅ (fu) = 0  → ∂f/∂t + (uᵢfₓ + uⱼf𝒚 + uₖf𝒛) = 0  
include("../../OZ_Tools/ExampleData.jl")
n = 10 # Number of spatial grid points
h = π/(n+1)
x = range(-π/2, step=h, length=n) 
y = range(-π/2, step=h, length=n) 
Δx = x[2]-x[1]
Δy = y[2]-y[1]
ord_deriv = 1
ord_approx = 2

# Δx = CenteredDifference{1}(ord_deriv, ord_approx, Δx, n)
# Δy = CenteredDifference{2}(ord_deriv, ord_approx, Δy, n)
Δuwx = UpwindDifference{1}(ord_deriv, ord_approx, Δx, n, 1)
Δuwy = UpwindDifference{2}(ord_deriv, ord_approx, Δy, n, 1)

u₀ = getGauss(n, n/7,  2)
heatmap(u₀, title="u₀")
##

# Uf(x,y) = -sin.(y) # Vector field U component
# Vf(x,y) =  sin.(x) .* sin.(y./2 .+ π/2 ) # Vector field V component
X = (x * ones(length(x))')';
Y = x * ones(length(x))';
# U = Uf(X,Y)
# V = Vf(X,Y)
U = 0
V = 1
##
function plotvector(x,y,u,v,subsample)
    l = 1:length(vec(X))
    lt = shuffle(l)
    ls = lt[1:subsample:end]
    quiver(vec(X)[ls],vec(Y)[ls],quiver=(vec(U)[ls],vec(V)[ls]))
end
# p3 = heatmap(x,x,U, title="U", c = :delta)
# p4 = heatmap(x,x,V, title="V", c = :delta)
# plotvector(X,Y,U,V,20)

bcx, bcy = RobinBC((1.0, 0.0, 0.0), (1.0, 0.0, 0.0), (Δx, Δy), 1, size(u₀)) # Defines Dirichlet BC with 0`
step(u,p,t) = (Δuwx*bcx*u).*U + (Δuwy*bcy*u).*V # Here is the definition of uₜ as an ODE
# step(u,p,t) = (Δuwx*bcx*u) + (Δuwy*bcy*u) # Here is the definition of uₜ as an ODE

t0 = 0.0
t1 = 2.0
Δt = 0.01
tsteps = Int(round(t1/Δt)) - 1
prob = ODEProblem(step, u₀, (t0, t1))
sol = solve(prob, saveat=Δt)
println("Done!")

##
for i=1:Int(round(tsteps/10)):tsteps
    t = sol.t[i]
    cur_true_sol = getGauss(n, n/7,  2, (-(V/Δx)*(Δt*i), -U*t*Δy))
    a = heatmap(cur_true_sol, c=:dense, title="Time t=$t true", clim=(-1,1))
    b = heatmap(sol[:,:,i]', c=:dense, title="Time t=$t UPW", clim=(-1,1))
    # plot(cur_true_sol[5,:], label="true", title = "Time t = $t")
    # p = plot!(sol[5,:,i], label="nwp")
    p = plot(a,b,layout = @layout[a b])
    xaxis!("X")
    yaxis!("Y")
    Plots.display(p)
    sleep(.1)
end
