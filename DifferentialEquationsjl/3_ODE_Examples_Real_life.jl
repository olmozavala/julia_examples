# Usings
using DifferentialEquations
using Plots
gr()

## =============  Population grow uₜ = k*u(1-u/N)- λ
# k → grow rate (Every generation how much the population increases)
# N → carrying capacity (The maximum 'allowed' population)
# λ → harvesting level (If we want to model harvesting in some way)
k = 1.1
N = 10^4
λ = 5.0

tspan = (0.0, 10.0)
u₀ = 10

# f(u,p,tspan) = k*u # Simply exponential grow
# f(u,p,tspan) = k*u*(1-u/N)  # Adding carrying capacity 
f(u,p,tspan) = k*u*(1-u/N) - λ # Adding carrying capacity and harvesting 
prob = ODEProblem(f, u₀, tspan)
sol = solve(prob) # Specify the error tolerance (compare and see)

plot(sol,lw=3,ls=:dash,title="Population Grow")

## =============  Predator-prey (Lotka-Volterra)
# xₜ = αx - βxy  # (Prey) α prey grow, β eating amount
# yₜ = δxy - γy  # (Predator) δ grow takes into account both species, γ desease rate

tspan = (0.0, 10.0)
u₀ = [20, 20] 

function lotka_volterra(du, u, p, t)
    α, β, δ, γ = p
    x, y  = u 
    du[1] = α*x - β*x*y
    du[2] = δ*x*y - γ*y
end

α = 4
β = .05
δ = .05
γ = 1
p = [α, β, δ, γ]

prob = ODEProblem(lotka_volterra, u₀, tspan, p)
sol = solve(prob) # Specify the error tolerance (compare and see)

plot(sol, title="Lotka Volterra", label=["1","2"])

## ============= Van der Pol oscillator uₜₜ = a(1-u²)uₜ - u
#https://en.wikipedia.org/wiki/Van_der_Pol_oscillator
tspan = (0.0, 10.0)

# 2*m*γ*u[1] - m*ω^2*u[2]
a = 1

function osci(du, u, p, t)
    # uₜₜ  = -u + a(1 - u²)uₜ
    # uₜ = z 
    # zₜ = -u + a(1 - u²)z
    # u → u[1] and uₜ = z = u[2]
    du[1] = u[2] # uₜ = z
    du[2] = a*(1 - u[1].^2)*u[2] - u[1]    # uₜₜ = -u + a(1 - u²)z
end

u₀ = [1, 1] # u₀ → 1 and uₜ → 1
prob = ODEProblem(osci, u₀, tspan)
sol = solve(prob) # Specify the error tolerance (compare and see)

l = @layout[a; b c]
a = plot(sol, title="Solution")
b = plot(sol.t, sol[1,:], title="u")
c = plot(sol.t, sol[2,:], title="u'", c=:Red)
plot(a, b, c, layout=l)
