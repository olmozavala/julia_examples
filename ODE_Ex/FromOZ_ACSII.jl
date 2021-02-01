# Usings
using DifferentialEquations
using Plots
gr()


## **************** Simplest case  du/dt = f(t,u) with f(t,u) = αu*************
println("Configuring problem...")
α = 1.1
f(u, p, t) = α*u
u₀ = 1/2
tspan = (0.0, 1.0)
typeof(tspan)
prob = ODEProblem(f, u₀, tspan)
# sol = solve(prob)
sol = solve(prob,reltol=1e-6,saveat=0.1) # Specify the error tolerance (compare and see)
# You need to run these lines manually or it wont show
plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
     xaxis="Time (t)",yaxis="u(t) (in μm)",label="ODE approx!") # legend=false
plot!(sol.t, t->u₀*exp(α*t),lw=3,ls=:dash,label="True Solution!")
print("Done!")

## ************* Systems of ODEs (lorenz) *************
function lorenz(t,u,du)
     du[1] = 10.0(u[2]-u[1])
     du[2] = u[1]*(28.0-u[3]) - u[2]
     du[3] = u[1]*u[2] - (8/3)*u[3]
end