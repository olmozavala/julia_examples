
using OrdinaryDiffEq
using ModelingToolkit
using DataDrivenDiffEq
using LinearAlgebra, DiffEqSensitivity, Optim
using DiffEqFlux, Flux
using Plots
gr()
using Statistics
# Set a random seed for reproduceable behaviour
using Random
Random.seed!(1234)

# Create a name for saving ( basically a prefix )
svname = "Scenario_1_"

## Data generation
function lotka!(du, u, p, t)
    α, β, γ, δ = p
    du[1] = α*u[1] - β*u[2]*u[1]
    du[2] = γ*u[1]*u[2]  - δ*u[2]
end

# Define the experimental parameter
tspan = (0.0f0,3.0f0)
u0 = Float32[0.44249296,4.6280594]
p_ = Float32[1.3, 0.9, 0.8, 1.8]
prob = ODEProblem(lotka!, u0,tspan, p_)
solution = solve(prob, Vern7(), abstol=1e-12, reltol=1e-12, saveat = 0.1)

# Ideal data
X = Array(solution)
t = solution.t

# Add noise in terms of the mean
x̄ = mean(X, dims = 2)
noise_magnitude = Float32(5e-2)
Xₙ = X .+ (noise_magnitude*x̄) .* randn(eltype(X), size(X))

plot(solution, alpha = 0.75, color = :black, label = ["True Data" nothing])
scatter!(t, transpose(Xₙ), color = :red, label = ["Noisy Data" nothing])
## ------------------ Define the UDE problem ---------------------------
# Gaussian RBF as activation
rbf(x) = exp.(-(x.^2))  # Why? 

# Define the network 2->5->5->5->2
U = FastChain(
    FastDense(2,5,rbf), FastDense(5,5, rbf), FastDense(5,5, rbf), FastDense(5,2)
)
# Get the initial parameters
p = initial_params(U)

# Define the hybrid model
function ude_dynamics!(du,u, p, t, p_true)
    û = U(u, p) # Network prediction
    du[1] = p_true[1]*u[1] + û[1]
    du[2] = -p_true[4]*u[2] + û[2]
end

# Closure with the known parameter
nn_dynamics!(du,u,p,t) = ude_dynamics!(du,u,p,t,p_)
# Define the problem
prob_nn = ODEProblem(nn_dynamics!,Xₙ[:, 1], tspan, p)

## Function to train the network
# Define a predictor (NODE) -------------------------
function predict(θ, X = Xₙ[:,1], T = t)
    Array(solve(prob_nn, Vern7(), u0 = X, p=θ,
                tspan = (T[1], T[end]), saveat = T,
                abstol=1e-6, reltol=1e-6,
                sensealg = ForwardDiffSensitivity()
                ))
end

# Simple L2 loss
function loss(θ)
    X̂ = predict(θ)
    sum(abs2, Xₙ .- X̂)
end

# Container to track the losses
# losses = Float32[]

# Callback to show the loss during training
callback(θ,l) = begin
    push!(losses, l)
    if length(losses)%50==0
        println("Current loss after $(length(losses)) iterations: $(losses[end])")
    end
    false
end

## Training ------------------------------

# First train with ADAM for better convergence -> move the parameters into a
# favourable starting positing for BFGS
res1 = DiffEqFlux.sciml_train(loss, p, ADAM(0.1f0), cb=callback, maxiters = 200)
println("Training loss after $(length(losses)) iterations: $(losses[end])")
# Train with BFGS
res2 = DiffEqFlux.sciml_train(loss, res1.minimizer, BFGS(initial_stepnorm=0.01f0), cb=callback, maxiters = 10000)
println("Final training loss after $(length(losses)) iterations: $(losses[end])")

# ===================================================================================
## Plot the losses
pl_losses = plot(1:200, losses[1:200], yaxis = :log10, xaxis = :log10, xlabel = "Iterations", ylabel = "Loss", label = "ADAM", color = :blue)
plot!(201:length(losses), losses[201:end], yaxis = :log10, xaxis = :log10, xlabel = "Iterations", ylabel = "Loss", label = "BFGS", color = :red)
# savefig(pl_losses, joinpath(pwd(), "plots", "$(svname)_losses.pdf"))
# Rename the best candidate
p_trained = res2.minimizer

## Analysis of the trained network
# Plot the data and the approximation
X̂ = predict(p_trained, Xₙ[:,1], t[1]:0.05f0:t[end])
# Trained on noisy data vs real solution
pl_trajectory = plot(t[1]:0.05f0:t[end], transpose(X̂), xlabel = "t", ylabel ="x(t), y(t)", color = :red, label = ["UDE Approximation" nothing])
scatter!(t, transpose(Xₙ), color = :black, label = ["Measurements" nothing])
# savefig(pl_trajectory, joinpath(pwd(), "plots", "$(svname)_trajectory_reconstruction.pdf"))

## ----------------
# n_tspan = (0.0f0,4.8f0)
# t = n_tspan[1]:0.1f0:n_tspan[end]
# prob_nn_oz = ODEProblem(nn_dynamics!, u0, n_tspan, res2.minimizer)
# sol_nn_oz = predict(p_trained, Xₙ[:,1], t)

# lot_oz = ODEProblem(lotka!, u0, n_tspan, p_)
# sol_oz = solve(lot_oz, saveat=0.1)

# plot(sol_oz)
# scatter!(t, transpose(sol_nn_oz))
##

# Ideal unknown interactions of the predictor
Ȳ = [-p_[2]*(X̂[1,:].*X̂[2,:])';p_[3]*(X̂[1,:].*X̂[2,:])']
# Neural network guess
Ŷ = U(X̂,p_trained)

pl_reconstruction = plot(t[1]:0.05f0:t[end], transpose(Ŷ), xlabel = "t", ylabel ="U(x,y)", color = :red, label = ["UDE Approximation" nothing])
plot!(t[1]:0.05f0:t[end], transpose(Ȳ), color = :black, label = ["True Interaction" nothing])
# savefig(pl_reconstruction, joinpath(pwd(), "plots", "$(svname)_missingterm_reconstruction.pdf"))

# Plot the error
pl_reconstruction_error = plot(t[1]:0.05f0:t[end], norm.(eachcol(Ȳ-Ŷ)), yaxis = :log, xlabel = "t", ylabel = "L2-Error", label = nothing, color = :red)
pl_missing = plot(pl_reconstruction, pl_reconstruction_error, layout = (2,1))
# savefig(pl_missing, joinpath(pwd(), "plots", "$(svname)_missingterm_reconstruction_and_error.pdf"))
pl_overall = plot(pl_trajectory, pl_missing)
# savefig(pl_overall, joinpath(pwd(), "plots", "$(svname)_reconstruction.pdf"))
## Symbolic regression via sparse regression ( SINDy based )

## Create a Basis
@variables u[1:2]
# Generate the basis functions, multivariate polynomials up to deg 5
# and sine
b = [polynomial_basis(u, 5); sin.(u)]
basis = Basis(b, u)
println(basis)

# Create an optimizer for the SINDy problem
opt = SR3(Float32(1e-2), Float32(0.1))
# Create the thresholds which should be used in the search process
λ = Float32.(exp10.(-7:0.1:5))
# Target function to choose the results from; x = L0 of coefficients and L2-Error of the model
g(x) = x[1] < 1 ? Inf : norm(x, 2)

# Test on ideal derivative data for unknown function ( not available )
println("SINDy on partial ideal, unavailable data")
Ψ = SINDy(X̂, Ȳ, basis, λ, opt, g = g, maxiter = 10000) # Succeed
println(Ψ)
print_equations(Ψ)

# Test on uode derivative data
println("SINDy on learned, partial, available data")
Ψ = SINDy(X̂, Ŷ, basis, λ,  opt, g = g, maxiter = 50000, normalize = true, denoise = true, convergence_error = Float32(1e-10)) # Succeed
println(Ψ)
print_equations(Ψ)

# Extract the parameter
p̂ = parameters(Ψ)
println("First parameter guess : $(p̂)")

# Just the equations
b = Basis((u, p, t)->Ψ(u, [1f0; 1f0], t), u)

# Retune for better parameters -> we could also use DiffEqFlux or other parameter estimation tools here.
cg = X̂
Ψf = SINDy(X̂, Ŷ, b, STRRidge(0.01f0), maxiter = 100, convergence_error = Float32(1e-18)) # Succeed
println(Ψf)
p̂ = parameters(Ψf)
println("Second parameter guess : $(abs.(p̂))")
println("True parameter : $(p_[2:3])")

# Define the recovered, hyrid model
function recovered_dynamics!(du,u, p, t, p_true)
    û = Ψf(u, p) # Network prediction
    du[1] = p_true[1]*u[1] + û[1]
    du[2] = -p_true[4]*u[2] + û[2]
end

# Closure with the known parameter
estimated_dynamics!(du,u,p,t) = recovered_dynamics!(du,u,p,t,p_)

estimation_prob = ODEProblem(estimated_dynamics!, u0, tspan, p̂)
estimate = solve(estimation_prob, Tsit5(), saveat = t)

# Plot
plot(solution)
plot!(estimate)

## Simulation

# Look at long term prediction
t_long = (0.0f0, 20.0f0)
estimation_prob = ODEProblem(estimated_dynamics!, u0, t_long, p̂)
estimate_long = solve(estimation_prob, Tsit5(), saveat = 0.1) # Using higher tolerances here results in exit of julia
plot(estimate_long)

true_prob = ODEProblem(lotka!, u0, t_long, p_)
true_solution_long = solve(true_prob, Tsit5(), saveat = estimate_long.t)
plot!(true_solution_long)

## Save the results
save(joinpath(pwd(), "results" ,"$(svname)recovery_$(noise_magnitude).jld2"),
    "solution", solution, "X", Xₙ, "t" , t, "neural_network" , U, "initial_parameters", p, "trained_parameters" , p_trained, # Training
    "losses", losses, "result", Ψf, "recovered_parameters", p̂, # Recovery
    "long_solution", true_solution_long, "long_estimate", estimate_long) # Estimation


## Post Processing and Plots

c1 = 3 # RGBA(174/255,192/255,201/255,1) # Maroon
c2 = :orange # RGBA(132/255,159/255,173/255,1) # Red
c3 = :blue # RGBA(255/255,90/255,0,1) # Orange
c4 = :purple # RGBA(153/255,50/255,204/255,1) # Purple

p1 = plot(t,abs.(Array(solution) .- estimate)' .+ eps(Float32),
          lw = 3, yaxis = :log, title = "Timeseries of UODE Error",
          color = [3 :orange], xlabel = "t",
          label = ["x(t)" "y(t)"],
          titlefont = "Helvetica", legendfont = "Helvetica",
          legend = :topright)

# Plot L₂
p2 = plot3d(X̂[1,:], X̂[2,:], Ŷ[2,:], lw = 3,
     title = "Neural Network Fit of U2(t)", color = c1,
     label = "Neural Network", xaxis = "x", yaxis="y",
     titlefont = "Helvetica", legendfont = "Helvetica",
     legend = :bottomright)
plot!(X̂[1,:], X̂[2,:], Ȳ[2,:], lw = 3, label = "True Missing Term", color=c2)

p3 = scatter(solution, color = [c1 c2], label = ["x data" "y data"],
             title = "Extrapolated Fit From Short Training Data",
             titlefont = "Helvetica", legendfont = "Helvetica",
             markersize = 5)

plot!(p3,true_solution_long, color = [c1 c2], linestyle = :dot, lw=5, label = ["True x(t)" "True y(t)"])
plot!(p3,estimate_long, color = [c3 c4], lw=1, label = ["Estimated x(t)" "Estimated y(t)"])
plot!(p3,[2.99,3.01],[0.0,10.0],lw=1,color=:black, label = nothing)
annotate!([(1.5,13,text("Training \nData", 10, :center, :top, :black, "Helvetica"))])
l = @layout [grid(1,2)
             grid(1,1)]
plot(p1,p2,p3,layout = l)

savefig(joinpath(pwd(),"plots","$(svname)full_plot.pdf"))
© 2021 GitHub, Inc.
Terms
Privacy
Security
Status
Docs
Contact GitHub
Pricing
API
Training
Blog
About
