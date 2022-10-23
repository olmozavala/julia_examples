## ===Code for advection Advection 1D  ∂f/∂t + ∇ ⋅ (fu) = 0  → ∂f/∂t + (uᵢfₓ + uⱼf𝒚 + uₖf𝒛) = 0   =================
using ModelingToolkit, DiffEqOperators, DifferentialEquations
using Plots
using Random
using DiffEqOperators, OrdinaryDiffEq
## 
gr()
include("../../OZ_Tools/ExampleData.jl")

## ---- Generating initial conditions u₀   
n = 100# Number of spatial grid points
Δx = 1.0
x = range(Δx, step=Δx, length=n) 
u₀ = getGaussNorm(n, n/8)
scatter(u₀, title="u₀")

## ---- Setup Finite difference on space (centered and upwind)
ord_deriv = 1
ord_approx = 2
Δ = CenteredDifference(ord_deriv, ord_approx, Δx, n) # Here we get the proper matrix for centered differencs
Δuw = UpwindDifference(ord_deriv, ord_approx, Δx, n, 1)
# display(Matrix(Δuw))
# println(size(Matrix(Δuw)))

## ---- Define boundary conditions
# Boundary conditions  au + b ∂u/∂n = c  
# bcx = RobinBC((1.0, 0.0, 0.0),(1.0, 0.0, 0.0), (Δx), 1, size(u₀)) # Direichlet = 0
# bcx = Neumann0BC(Δx)
# bcx = Dirichlet0BC(Float64)
bcx = PeriodicBC(Float64)  # x[-1] = x[end] and x[end+1] = x[0]

scatter(x, u₀, label="Scalar field")
# Plot the boundary conditions with u₀
display(scatter!(transpose([x[1]-Δx x... x[end]+Δx]), bcx*u₀, label="u₀ and boundary condition", 
                title="u₀ and boundary condition",
                markersize=2, markerstrokewidth=0, mcolor=:red, c=:red))
# Example of 'manual' BC
# display(scatter!(transpose([x[1]-Δx x... x[end]+Δx]), [0 u₀... 0]', label="Manual BC", markersize=1, markerstrokewidth=0, mcolor=:green, c=:red))
## -------- Define vector field (in this case it is static, but it could vary through time)
v = 1 # lets say speed in m per sec
U = v*ones(n) # Vector field moving homogeneously to one side
# U = -cos.(π*(x/x.len))# Vector field moving to both sides from the middle
plot(x, U, c=:green, label="Vector field", legend=:bottomright)


## -------- Defining the ODE
# === Advection Advection 1D  ∂f/∂t + ∇ ⋅ (fu) = 0  → ∂f/∂t + uᵢfₓ = 0 → ∂f/∂t = -uᵢfₓ  =================
# Define as a function with respect to U
step(f,p,t) = -1*Δ*bcx*(f.*U)  # Define the ODE with CenteredDifference
stepuw(f,p,t) = -1*Δuw*bcx*(f.*U)  # Define the ODE with UpwindDifference

## Solving equation
t0 = 0.0
t1 = 100
Δt = 1 # Lets say seconds
tsteps = Int(round(t1/Δt)) - 1
println("Solving with Δt=$Δt and Δx=$Δx")
prob = ODEProblem(step, u₀, (t0, t1))
sol = solve(prob, saveat=Δt);
probuw = ODEProblem(stepuw, u₀, (t0, t1))
soluw = solve(probuw, saveat=Δt);
println("Done!")

## Same thing but with an integrator
# integrator = init(probuw, Tsit5(), dt=Δt)
# # for i in 1:tsteps
# for i in 1:5
#     plot(x, -sin.(x .+ v*(Δt*(i-1))), title="Integrator", label="True t=$(Δt*(i-1))", c=:green, lw=4)
#     display(plot!(x, integrator.u, c=:yellow, mstyle=:dash, label="nm t=$(integrator.t))"))
#     step!(integrator, Δt, true)
#     sleep(1)
# end

## Plotting the result
println("Plotting...")
plot(u₀, title="Solution t=0")
lw = 4
x = range(Δx, step=Δx, length=n) 
# for i=1:Int(round(tsteps/10)):tsteps
for i=1:2:tsteps
    # p = plot(-sin.(x .+ v*(Δt*(i-1))), title="Solution True t=$(soluw.t[i])", label="True", c=:green, lw=lw)
    u_true = circshift(getGaussNorm(n, n/8, 0.0, 1) ,(v/Δx)*(Δt*(i-1)))
    # println((v/Δx)*(Δt*(i-1)))
    p = plot(x, u_true, title="Solution True t=$(soluw.t[i]) Speed: $v m/s", label="True", c=:green, lw=lw)
    plot!(sol[:,i], c=:blue, label="CD")
    plot!(soluw[:,i], c=:red, label="UW", linestyle=:dash, lw=lw)
    Plots.display(p)
    sleep(.1)
end
println("Done")


## ======================= 2D ================================
## ======================= 2D ================================
## ========= Advection 2D  ∂f/∂t + ∇ ⋅ (fu) = 0  → ∂f/∂t + (uᵢfₓ + uⱼf𝒚) = 0  
# Define initial condition of f(x,y)
include("../../OZ_Tools/ExampleData.jl")
n = 100 # Number of spatial grid points
Δx = 1.0
Δy = 1.0
x = range(0, step=Δx, length=n) 
y = range(0, step=Δy, length=n) 
ord_deriv = 1
ord_approx = 2

# Lo que esta entre llaves es el axis
# Δx = CenteredDifference{1}(ord_deriv, ord_approx, Δx, n)
# Δy = CenteredDifference{2}(ord_deriv, ord_approx, Δy, n)
Δuwx = UpwindDifference{1}(ord_deriv, ord_approx, Δx, n, 1)
Δuwy = UpwindDifference{2}(ord_deriv, ord_approx, Δy, n, 1)

u₀ = getGaussNorm(n, n/30,  0.0, 2)
heatmap(u₀, title="u₀")

## Define velocity vector filed U,V
# ---------- Not sure
# Uf(x,y) = -sin.(y) # Vector field U component
# Vf(x,y) =  sin.(x) .* sin.(y./2 .+ π/2 ) # Vector field V component
# U = Uf(X,Y)
# V = Vf(X,Y)
# ---------- Moving only horizontally
X = (x * ones(length(x))')';
Y = y * ones(length(y))';
U = ones(size(X))*1  # U=1 everywhere (moving to the right)
V = zeros(size(X))   # V=0 everywhere (not moving)

# Plotting the velocity vector
subsample = 10 
function plotvector(x,y,u,v,subsample)
    l = 1:length(vec(X))
    lt = shuffle(l)
    ls = lt[1:subsample:end]
    quiver(vec(X)[ls],vec(Y)[ls],quiver=(vec(U)[ls],vec(V)[ls]))
end
p3 = heatmap(x,x,U, title="U", c = :delta)
p4 = heatmap(x,x,V, title="V", c = :delta)
plot(p3,p4, layout=@layout [a b])
# plotvector(X,Y,U,V,20)

##
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

## Plotting the results
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
