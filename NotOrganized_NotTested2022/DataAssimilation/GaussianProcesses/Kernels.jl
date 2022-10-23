
using Distributions
using Plots
using GaussianProcesses
using Random
using Optim
using Distances
gr()

# Squared Exponential covariance function
k₁(x,y) = e.^(-(1/2)*evaluate(Euclidean(), x, y())
kern = SE(0.0,0.0)

