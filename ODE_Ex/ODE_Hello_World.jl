# Usings
using DifferentialEquations
using Plots
gr()

## **************** Simplest case  du/dt = f(t,u) with f(t,u) = αu*************
println("Configuring problem...")
α = 1.1
Δx = 0.1
f(u, p, t) = α*u
u₀ = 1/2
tspan = (0.0, 1.0)
x = tspan[1]:Δx:tspan[end]
typeof(tspan)
prob = ODEProblem(f, u₀, tspan)
# sol = solve(prob)
sol = solve(prob,reltol=1e-6,saveat=Δx) # Specify the error tolerance (compare and see)
# You need to run these lines manually or it wont show
plot(x, t->u₀*exp(α*t), lw=5, label="True Solution!")
plot!(sol,lw=3,ls=:dash,title="Solution to the linear ODE with a thick line",
     xaxis="Time (t)",yaxis="u(t) (in μm)",label="ODE approx!") # legend=false
print("Done!")

## ******* Same solved with finite differences *******
using DiffEqOperators, LinearAlgebra
order = 2
deriv = 1
Δx = 0.1
N = length(x)

# Setting  the RHS
b = zeros(N)
b[1] = u₀

# Manually building A with BE
A = zeros(N, N) 
A[1] = 1
for i = 2:N
     A[i,i-1] = -α*Δx - 1
     A[i,i] = 1
end
display(A)

ufd = A\b
plot(x, t->u₀*exp(α*t), lw=5, label="True Solution!")
scatter!(x, ufd, label="FD Manually 1 order", legend=:topleft)

## Manually building A with central differences
A = zeros(N, N) 
A[1,1] = 1
A[2,2-1] = -α*Δx - 1
A[2,2] = 1
for i = 3:N
     A[i,i-2] = -1
     A[i,i-1] = -α*2*Δx
     A[i,i] = 1
end

display(A)
ufd = A\b
plot(x, t->u₀*exp(α*t), lw=5, label="True Solution!")
scatter!(sol.t, ufd, label="FD Manually 2nd order", legend=:topleft)