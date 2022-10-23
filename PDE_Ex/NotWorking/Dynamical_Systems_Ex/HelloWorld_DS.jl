# Docs https://juliadynamics.github.io/DynamicalSystems.jl/latest/
# Tutorials https://github.com/JuliaDynamics/JuliaDynamics/tree/master/tutorials/Youtube_JuliaLang_tutorial
using DynamicalSystems
using Plots
# plotly()
gr()

# =============== First examples only DiscreteDS and ContinuousDS ====================
# Example 1 Henon map (discreate DS):   https://en.wikipedia.org/wiki/H%C3%A9non_map
#       xₙ₊₁ = 1 - ax²ₙ + yₙ
#       yₙ₊₁ = bxₙ 

# h_eom(x, p, t) = SVector{2}(1 - p[1]x[1]^2 + x[2], p[2]*x[1])
# or
function h_eom(du, x, p, t)
    du[1] = 1 - p[1]x[1]^2 + x[2]
    du[2] = p[2]*x[1]
    return nothing
end
state = zeros(2)  # Initial state
p = [1.4, 0.3]   # Parameters
tsteps = 500

henon = DiscreteDynamicalSystem(h_eom, state, p)
display(henon) # Display what was created
# tr = trajectory(henon, tsteps)
tr = trajectory(henon, tsteps, 0.01rand(2)) # Different initial state

scatter(tr[:,1], tr[:,2], title="Example 1 Henon map")

## Example 1 Henon-Heiles (continuous DS):   https://en.wikipedia.org/wiki/H%C3%A9non%E2%80%93Heiles_system
function hheom!(du, u, p, t)
    du[1] = u[3]
    du[2] = u[4]
    du[3] = -u[1] - 2u[1]*u[2]
    du[4] = -u[2] - (u[1]^2 - u[2]^2)
    return nothing
end
state = 0.1rand(4)  # Initial state
state = [0, -0.25, 0.42081, 0]
tsteps = 100

hh= ContinuousDynamicalSystem(hheom!, state, nothing)
display(hh) # Display what was created
trhh = trajectory(hh, tsteps)

plot(trhh[:,1], trhh[:,2], title="Example 2 HH")

## =============== Orbit diagram  ========
# Analyze the Long term behaviour of some of the variables

# Exmaple with the logistic map xₙ₊₁ = rxₙ(1-xₙ)  https://en.wikipedia.org/wiki/Logistic_map
include("PlotToolsDS.jl")
logimap = Systems.logistic() # Define DS problem
i = 1 # Which variable index we want to show
p_index = 1 # Which parameter we want to analyze 
pvalues = 3.4:0.001:3.75  # parameter values
n = 200 # how many values to save for each 'run'
Ttr = 1 # how many iterations to 'skip'. Not clear what is this??????????
output = orbitdiagram(logimap, i, p_index, pvalues; n = n, Ttr = Ttr)
typeof(output)
plot_od(output, pvalues, n)

## =============== Poincare Surface of Section ========
# Record hyperplanes for desired areas https://en.wikipedia.org/wiki/Poincar%C3%A9_map
lor = Systems.lorenz()
tr = trajectory(lor, 100.0, dt = 0.005, Ttr = 50.0)
x, y, z = columns(tr)
plot(x,y,z, leg=false, title="Lorenz attractor", 
        html_output_format=:png, size=(1000,1000))

## Define and plot the plane
plane = (2, 0.0) # when 2nd variable crosses 0.0
psos_chaotic = poincaresos(lor, plane, 2000.0, Ttr = 100.0)

scatter(psos_chaotic[:, 1], psos_chaotic[:, 3],
        markersize=0.15, markeralpha = 0.15, markercolor=:black,
        leg=false, title="Lorenz attractor Poincare section", 
        html_output_format=:png, size=(500,500))


## =============== Lyapunov exponents ==============
using LinearAlgebra: norm
henon = Systems.henon()
# Define two starting points that are close enough
tr1 = trajectory(henon, 100)
u2 = get_state(henon) + (1e-9 * ones(dimension(henon))) # Close initial state
tr2 = trajectory(henon, 100, u2)

λ = lyapunov(henon, 5000) # Compute Lyapunov second argument is time to evolve
λs = lyapunovs(henon, 5000) # Compute all Lyapunovs

## Plot results
p1 = plot(tr1[:, 1], alpha = 0.5, title="Diverging trajectories")
plot!(tr2[:, 1], alpha = 0.5)

d = [norm(tr1[i] - tr2[i]) for i in 1:length(tr2)]
p2 = plot(d, yaxis=:log, title="Their distance")
plot!(collect(0:50), d[1] .* exp.(collect(0:50) .* λ))

plot(p1,p2,layout=(2,1),legend=false)