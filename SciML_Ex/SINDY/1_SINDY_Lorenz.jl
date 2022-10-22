using DataDrivenDiffEq
using ModelingToolkit
using OrdinaryDiffEq

using LinearAlgebra
using Plots
gr()
##

# Create a test prob
function lorenz(u,p,t)
    x, y, z = u
    ẋ = 10.0*(y - x)
    ẏ = x*(28.0-z) - y
    ż = x*y - (8/3)*z
    return [ẋ, ẏ, ż]
end

u0 = [-8.0; 7.0; 27.0]
p = [10.0; -10.0; 28.0; -1.0; -1.0; 1.0; -8/3]
tspan = (0.0,100.0)
dt = 0.001
prob = ODEProblem(lorenz,u0,tspan)
sol = solve(prob, Tsit5(), saveat = dt, atol = 1e-7, rtol = 1e-8)
x = [x[1] for x in sol.u][1:50:50000]
y = [x[2] for x in sol.u][1:50:50000]
z = [x[3] for x in sol.u][1:50:50000]
## ---- Plot Lorenz -------
plt = plot3d(
    1,
    xlim = (-30, 30),
    ylim = (-30, 30),
    zlim = (0, 60),
    title = "Lorenz Attractor",
    marker = 2,
)
# build an animated gif by pushing new points to the plot, saving every 10th frame
@gif for i=1:size(x)[1]
    push!(plt, x[i], y[i], z[i])
end every 10
##
# -------------- Define problem (from your data) ----------------------
X = sol[:,:]
ts = sol.t
prob = ContinuousDataDrivenProblem(X, ts)

## -------------- Generate Basis ----------------------
@variables u[1:2]
# h = Num[sin.(w[1].*u[1]);cos.(w[2].*u[1]); polynomial_basis(u, 5); c]
h = polynomial_basis(u, 4)
basis = Basis(h, u)

## --------- --------
sampler = DataSampler(Batcher(n = 50, shuffle = true, repeated = true))
λs = exp10.(-10:0.1:-1)
opt = STLSQ(λs)
res = solve(prob, basis, opt, progress = true, sampler = sampler, 
                    denoise = false, normalize = true, maxiter = 5000)
print("Done!")
print(res)
# println(result(res))