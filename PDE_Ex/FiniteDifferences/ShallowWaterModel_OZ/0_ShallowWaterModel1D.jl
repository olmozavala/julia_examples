using ModelingToolkit, DiffEqOperators, DifferentialEquations
using Plots
using Random
using DiffEqOperators, OrdinaryDiffEq
gr()
include("/home/olmozavala/Dropbox/TutorialsByMe/Julia/Julia_Examples/OZ_Tools/ExampleData.jl")


## ========= Simple ShallowWaterModel 1D ================
# channel_size = 150  # Kilometers
# Δx = 1000.0  # Meters
# n = Int(round(channel_size*1000/Δx)) # Number of spatial grid points
n = 50
Δx = 2*π/n
x = range(Δx, step=Δx, length=n) 
waveheight = 20
# h₀ = waveheight*getGaussShape(n, n/8, 1, 0.0) 
# h₀ = zeros(n)
h₀ = waveheight*sin.(x)

# Boundary conditions  au + b ∂u/∂n = c  
# bc = RobinBC((1.0, 0.0, 0.0),(1.0, 0.0, 0.0), (Δx), 1, size(u₀)) # Direichlet = 0
bcn = Neumann0BC(Δx)
bcd = Dirichlet0BC(Float64)
bcp = PeriodicBC(Float64)

scatter(x, h₀, label="SSH",  markersize=9, markerstrokewidth=0)
scatter!(transpose([x[1]-Δx x... x[end]+Δx]), bcp*h₀, label="BC Periodic", markersize=6, markerstrokewidth=0, mcolor=:red, c=:red)
scatter!(transpose([x[1]-Δx x... x[end]+Δx]), bcd*h₀, label="BC Dirichlet0", markersize=4, markerstrokewidth=0, mcolor=:black, c=:black)
p = scatter!(transpose([x[1]-Δx x... x[end]+Δx]), bcn*h₀, label="BC Neumann0", markersize=2, markerstrokewidth=0, mcolor=:cyan, c=:black)
xaxis!("meters")
yaxis!("meters")
display(p)

## We set f in the following way [u h]
ord_deriv = 1
ord_approx = 2
Δ = CenteredDifference(ord_deriv, ord_approx, Δx, n)

# We first set u and then SSH
f₀ = [zeros(n)...  h₀...]
bcu = bcp # Set u BC to zero
bcssh = bcp  # Set ssh BC to periodic
function step(du,f,p,t)
    n = p[1]
    u = f[1:n]
    ssh = f[n+1:2*n]

    # if t < 151
    #     periodo = 150.0 # Metros
    #     # bcssh = DirichletBC(sin(t*2*π/periodo), 0.0)
    #     bcssh = DirichletBC(1.0, 0.0)
    # else
    #     # After t = 151 we add periodic BC
    #     # bcssh = bcp
    #     bcssh = DirichletBC(-1.0, 0.0)
    # end
    g = 9.81 # Gravity mt/s²
    n = p[1]
    H = 1 # Depth of the channel
    du[1:n] = -g*Δ*bcssh*ssh # Solving for u
    du[n+1:2*n] = -H*Δ*bcu*u # Solving for h
end

# Solving equation
t0 = 0.0
# t1 = 3600*8.0 #  
t1 = 2.0
Δt =  .02 # Seconds

# CFL test don't pay attention to it
println("CFL condition for v = 1. Δt should be ≤ $(sum(1/(ones(n).*(1/Δx))))")

tsteps = Int(round(t1/Δt)) - 1
println("Solving with Δt=$Δt and Δx=$Δx")
p = (n) # Parameters (just the size of n)
prob = ODEProblem(step, f₀, (t0, t1), p)
sol = solve(prob, saveat=Δt);
println("Done!")

##  
println("Plotting...")
lw = 4
for i=1:Int(round(tsteps/24)):tsteps
# for i=1:100:length(sol.t)
# for i=1:10:Int(t1)
    a = plot(sol.u[i][1:n], c=:blue, label="U", title="Time = $(sol.t[i]) sec", ylim=[-4*waveheight,4*waveheight])
    b = plot(sol.u[i][n+1:2*n], c=:green, label="SSH", title="Time = $(sol.t[i]) sec", ylim=[-waveheight, waveheight])
    p = plot(a,b, layout=@layout[a b], size=(1200,600))
    xaxis!("Meters")
    yaxis!("Meters")
    Plots.display(p)
    sleep(.2)
end

println("Done")