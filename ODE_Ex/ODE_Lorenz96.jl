
# Usings
using DifferentialEquations
using Plots
gr()


## ************* Lorenz N=3 *************
function lorenz(du, u, p, t)
    du[1] = p[1]*(u[2]-u[1])
    du[2] = u[1]*(p[2]-u[3]) - u[2]
    du[3] = u[1]*u[2] - p[3]*u[3]
end

u₀ = [1.0;0.0;0.0]
tspan = (0.0,100.0)
p = [10.0, 28.0, 8/3]
prob = ODEProblem(lorenz,u₀,tspan, p)
sol = solve(prob)

plot(sol,vars=(1,2,3))


## ************* Lorenz N=40*************
N = 40
F = 8

function lorenz96(du, u, p, t)
    F = p[1] 
    N = p[2]
    du[1] = (u[2] - u[N-1])*u[N] - u[1] + F
    du[2] = (u[3] - u[N])  *u[1] - u[2] + F
    du[N] = (u[1] - u[N-2])*u[N-1] - u[N] + F
    for i=3:N-1
        du[i] = (u[i+1] - u[i-2])*u[i-1] - u[i] + F
    end
end

u₀ = F * ones(N)
u₀[1] += 0.01
tspan = (0.0, 30.0)
tsteps = range(tspan[1], tspan[2], length = Int(ceil(tspan[2]/.05)))
p = [F, N]
prob = ODEProblem(lorenz96,u₀,tspan, p)
sol = solve(prob, saveat=tsteps)

# plot(sol,vars=(1,2,3))
heatmap(sol[:,:])

## Compute Lyapunov exponents
u2₀ = u₀ += .001 * ones(N)
prob2 = ODEProblem(lorenz96,u2₀,tspan, p)
sol2 = solve(prob2, saveat=tsteps)

l = @layout[a;b]
p1 = plot(sol[1,1:100])
plot!(sol2[1,1:100])

until_t = 100
start = 10
p3 = plot([norm(sol[1,1:i] - sol2[1,1:i]) for i=start:until_t])
for j=2:N
    plot!([norm(sol[j,1:i] - sol2[j,1:i]) for i=start:until_t])
end

plot(p1, p3, layout = l, legend=false)


# Compute maximal Lyapunov
t = 600 # last time
λ = 0
for i=1:t
    norma = norm(sol[:,1] - sol2[:,1])
    normt = norm(sol[:,i] - sol2[:,i])
    λ += log(normt/norma)
end
λ *= 1/N