using Distributions
using Plots
using GaussianProcesses
using Random
using Optim
gr()

# https://juliahub.com/docs/GaussianProcesses/izlaf/0.12.0/Regression/
# ========= Simulate Training data 
n=10;                          #number of training points
x = 2π * rand(n);              #predictors
y = sin.(x) + 0.05*randn(n);

scatter(x,y, title="Simulated data")

## ======= Make your GaussianProcess using a kernel function and a mean value
mZero = MeanZero() #Zero mean function
#Squared exponential kernel (l,σ) (note that hyperparameters are on the log scale)
# l → length scale
# σ → signal standar deviation
kern = SE(0.0,0.0)
logObsNoise = -1.0 # log standard deviation of observation noise (this is optional)
gp = GP(x,y,mZero,kern,logObsNoise)       #F

# ======  After fitting to a GP, we can predict new values of unobserved points
μ, σ² = predict_y(gp,range(0,stop=2π,length=100));
plot(gp; xlabel="x", ylabel="y", title="Gaussian process", legend=false)      # Plot the GP

## ===== Optimizing parameters 
res = optimize!(gp) # Compute new parameters
display(res.minimizer)
# gp has already been modified so we can just plot the new results
a = plot(gp; xlabel="x", ylabel="y", title="Gaussian process", legend=false)      # Plot the GP
## If we want to 'probe' the results we redo the problem
r = res.minimizer 
println("Optimized parameters: \n observation noise = $(r[1]) ")
println("Length scale = $(r[2]) \n signal standard deviation = $(r[3]) \n ")
kern = SE(r[2], r[3])
logObsNoise = r[1] 
ngp = GP(x,y,mZero,kern,logObsNoise)       #F
μ, σ² = predict_y(ngp,range(0,stop=2π,length=100));
l = @layout[a b]
b = plot(ngp; xlabel="x", ylabel="y", title="Gaussian process", legend=false)      # Plot the GP
plot(a, b, layout=l)




