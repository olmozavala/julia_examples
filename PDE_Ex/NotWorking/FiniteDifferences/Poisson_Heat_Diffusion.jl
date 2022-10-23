using ModelingToolkit, DiffEqOperators, DifferentialEquations
using Plots
using DiffEqOperators, OrdinaryDiffEq
gr()
# # Heat/Diffusion Equation 1D uₜ = Δu  or uₜ(t,x,y) = uₓₓ + u𝐲𝐲 # Heat
#
# This example demonstrates how to combine `OrdinaryDiffEq` with `DiffEqOperators` to solve a time-dependent PDE.
# We consider the heat equation on the unit interval, with Dirichlet boundary conditions:
# ∂ₜu = Δu
# u(x=0,t)  = a
# u(x=1,t)  = b
# u(x, t=0) = u₀(x)
#
# For `a = b = 0` and `u₀(x) = sin(2πx)` a solution is given by:
u_analytic(x, t) = sin(2*π*x) * exp(-t*(2*π)^2)
#
# We want to reproduce it numerically


## ======= Define grid and boundary conditions ===========
nknots = 100 # Number of spatial grid points
h = 1.0/(nknots+1)
x = range(h, step=h, length=nknots) # This is in fact our grid points in space
ord_deriv = 2
ord_approx = 2

Δ² = CenteredDifference(ord_deriv, ord_approx, h, nknots)
bc = Dirichlet0BC(Float64) # Defines Dirichlet BC with 0
display(Matrix(Δ²))  # Displays the BC

t0 = 0.0
t1 = 0.03
u0 = sin.(2*π*x)

scatter(x, u0, ms=5)
s1 = scatter!([0 x... x[end]+h]', bc*u0, title="u₀ with 0 Dirichlet BC", c=:red, ms=2) # Here we plot u₀ with the boundary conditions

## ===== Solve heat equation as an ODE on the right hand side 
println("Solving equation as an ODE...")
step(u,p,t) = Δ²*bc*u  # Here is the definition of uₜ as an ODE (Diffusion)
# Another example with time dependent boudnary conditions
#function step(u,p,t)
#    bc = DirichletBC(cos(t*1.5), sin(t))
#    Δ*bc*u
#end
prob = ODEProblem(step, u0, (t0, t1))
# alg = KenCarp4()
sol = solve(prob, saveat=.001)
println("Done!")

l = @layout[a b]
h1 = heatmap(hcat(sol.u...))
xaxis!("Time")
yaxis!("Space")
plot(s1, h1, layout=l)
## ===== Same idea but for a diffusion equation (same but with k)
κ = 20
step(u,p,t) = κ*Δ*bc*u  # Here is the definition of uₜ as an ODE
# Another example with time dependent boudnary conditions
#function step(u,p,t)
#    bc = DirichletBC(cos(t*1.5), sin(t))
#    Δ*bc*u
#end
prob = ODEProblem(step, u0, (t0, t1))
alg = KenCarp4()
sol = solve(prob, saveat=.001)

l = @layout[a b]
h1 = heatmap(sol[:,:]', title="Diffusion Solution k = $κ")
xaxis!("Space")
yaxis!("Time")
plot(s1, h1, layout=l)


## ========= Heat Equation 2D
nknots = 20 # Number of spatial grid points
h = 1.0/(nknots+1)
x = range(h, step=h, length=nknots) # This is in fact our grid points in space
y = range(h, step=h, length=nknots) # This is in fact our grid points in space
ord_deriv = 2
ord_approx = 2

Δxx = CenteredDifference{1}(ord_deriv, ord_approx, h, nknots)
Δyy = CenteredDifference{2}(ord_deriv, ord_approx, h, nknots)
Δ = Δxx + Δyy
# bcx = Dirichlet0BC{1}() # Defines Dirichlet BC with 0
# bcy = Dirichlet0BC{2}() # Defines Dirichlet BC with 0
# bc = Dirichlet0BC(Float64) # Defines Dirichlet BC with 0
Qx, Qy = Neumann0BC(Float64, (h, h), 4, (length(x), length(y)))
Q = compose(Qx, Qy)

t0 = 0.0
t1 = 0.03
f(x,y) = sin.(2*π*x) + sin.(2*π*y)
u₀ = zeros(length(x), length(y))
for i = 1:length(x), j = 1:length(y)
    u₀[i,j] = f(x[i], y[j])
end
heatmap(u₀)

step(u,p,t) = Δ*Q*u # Here is the definition of uₜ as an ODE
# Second example with time dependent boudnary conditions
#function step(u,p,t)
#    #bc = DirichletBC(cos(t*1.5), sin(t))
#    bc = Dirichlet0BC(Float64)
#    Δxx*bc*u + Δxy*bc*u # Here is the definition of uₜ as an ODE
#end
#prob = ODEProblem(step, u₀, (t0, t1))
prob = ODEProblem(step, u₀, (t0, t1))
sol = solve(prob, saveat=.001)

for i=1:size(sol)[1]
    p = heatmap(sol[i,:,:]', title="Solution t=$i")
    xaxis!("X")
    yaxis!("Y")
    Plots.display(p)
end
