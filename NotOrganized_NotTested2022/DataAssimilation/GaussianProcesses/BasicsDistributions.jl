using Distributions
using Plots
gr()

d = Normal()
scatter(rand(d,100)