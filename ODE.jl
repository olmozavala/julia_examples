# Usings
using DifferentialEquations
using Plots
gr()

# **************** Simplest case *************
## =================== Numerically solving ODE ====
println("Configuring problem...")

a = 0.98
f(u, p, t) = a*u
u₀ = 1.0
tspan = (-1.0, 1.0)
typeof(tspan)
prob = ODEProblem(f, u₀, tspan)
sol = solve(prob)
println("Done!")

## =================== Plotting/Printing results

println("Making plot...")
println("Solution values: $sol")
plot(Array(sol.t), real)
plot(sol.t, sol.u)
# real = exp.(a .* Array(sol.t))
# print(real)
# print(Array(sol.t))
# println("Done!")