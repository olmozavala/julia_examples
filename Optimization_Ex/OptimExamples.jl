using Optim
using Plots
pyplot()
## https://julianlsolvers.github.io/Optim.jl/stable/#user/minimization/

# === Hello world rosenbrock function
rosenbrock(x) =  (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2

x = -1.5:.02:1.5
n = length(x)
Z = zeros(n,n)
for i=1:n, j=1:n
    Z[i,j] = rosenbrock([x[i],x[j]])
end
# heatmap(x, x, Z, c=:rainbow)
plot(x, x, Z, c=:rainbow, st=:surface, camera=(30,50))

## Plot results
result = optimize(rosenbrock, [1.0 -1.0], BFGS())
# List of functions here:
# https://julianlsolvers.github.io/Optim.jl/stable/#user/minimization/#complete-list-of-functions
scatter!([result.minimizer[1]], [result.minimizer[2]], [0.0], c=:red, ms=10)