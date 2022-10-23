using LinearAlgebra
using Plots
using Random, Distributions
using StatsBase
gr()
# Base.run(`clear`)

## ============ Covariance Matrix ==========
## Generate some data
tsteps = 10000
# Standar deviations
σ₀ = 1   
σ₁ = .5 
T₀ = rand(Normal(0,σ₀), tsteps)
T₁ = rand(Normal(0,σ₁), tsteps)

scatter(T₀, T₁, label="Obs", msw = 0, mc = :blue, ms=2)

## Compute covariance matrix
m₀ = mean(T₀)
m₁ = mean(T₁)
# m, c = mean_and_cov([T₀ T₁])

# -- Manually cov matrix ---
A = zeros(2,2)
A[1,1] = ((T₀ .- m₀)'*(T₀ .- m₀))/tsteps
A[1,2] = ((T₁ .- m₁)'*(T₀ .- m₀))/tsteps
A[2,1] = ((T₀ .- m₀)'*(T₁ .- m₁))/tsteps
A[2,2] = ((T₁ .- m₁)'*(T₁ .- m₁))/tsteps

# -- Automatic cov matrix ---
c = cov([T₀ T₁])
# display(c)
# display(A)

## ================ Distance correlation functions =======
# From (Page 310) https://www.hycom.org/attachments/428_cummings-chap13.pdf

# Just making a function to plot some points of the resulting matrices
function plotResults(SOAR, Gauss, w, h)
    lay = @layout[a b]
    a = heatmap(SOAR, title="SOAR", aspect_ratio=:equal)
    b = heatmap(Gauss, title="Gauss", aspect_ratio=:equal)
    display(plot(a,b, layout=lay, size=(1000,500)))

    for i=Int(round(w/2)):Int(round(w * round(h/4))):w*h # Plot four values at the middle
    # for i=1:w*10:w*h
    # for i=1:2:w*h
        lay = @layout[a b]
        cur_SOAR = reshape(SOAR[i,:], (w,h))
        cur_G = reshape(Gauss[i,:], (w,h))
        a = heatmap(transpose(cur_SOAR), title="SOAR", aspect_ratio=:equal)
        b = heatmap(transpose(cur_G), title="Gauss", aspect_ratio=:equal)
        display(plot(a,b, layout=lay, size=(600,600)))
        println("Plotting $i")
    end
end

## Second order auto regressive function SOAR
l = 10 # Number of locations per axis (lat,lon)
n = l*l # Total number of correlation values
SOAR = zeros(n, n) # Correlation matrix with size l² 
Gauss = zeros(n, n) # Correlation matrix with size l² 

# Here we build SOAR row by row (each row has all the spatial grid points)
for i=1:l, j=1:l, row=1:n 
    col = (i-1)*l + j
    a = [ceil(row/l), (row-1)%l+1] # First location we are comparing
    b = [i, j] # second location
    # println("$row x $col  --> $a,$b")
    # ---------- SOAR ----
    s = norm(a-b) 
    SOAR[row, col] = (1 + s)*exp(-s)
    # ---------- Gaussian ----
    Gauss[row, col] = exp(-s)
end
plotResults(SOAR, Gauss, l, l)

## Custom size correlation function with proper Rossby deformation taked into account.
using Distances

LON = [-10,10]
LAT = [0,90]
res = 1.0 # Resolution on degrees
nₕ = Int(ceil((maximum(LON)-minimum(LON))/res)) + 1
nᵥ = Int(ceil((maximum(LAT)-minimum(LAT))/res)) + 1
println("Size of matrix: $nₕ x $nᵥ")

# Second order auto regressive function SOAR
n = nₕ*nᵥ # Total number of correlation values
SOAR = zeros(n, n) # Correlation matrix with size l² 
Gauss = zeros(n, n) # Correlation matrix with size l² 

# Here we build SOAR row by row (each row has all the spatial grid points)
sop = 1
max_angle = 30
println("Building matrix...")
for row=1:n, c_lat=1:nᵥ, c_lon=1:nₕ
    col = (c_lat-1)*nₕ + c_lon 
    a = [floor(row/nₕ),  row-(floor(row/nₕ)*nₕ)-1] # First location we are comparing
    b = [c_lat-1, c_lon-1] # second location indexes

    p1 = [LON[1]+a[2]*res, LAT[1]+a[1]*res]
    p2 = [LON[1]+b[2]*res, LAT[1]+b[1]*res]
    
    # println("$p1 -- $p2")
    if abs(p1[2] - p2[2]) < max_angle  && abs(p1[1] - p2[1]) < max_angle
        s = haversine(p1, p2, 6372.8)/500
        # println("$p1, $p2 --> $s")
        # ---------- SOAR ----
        SOAR[row, col] = (1 + s)*exp(-s)
        # SOAR[row, col] = (1 + s/L)*exp(-s/L)
        # ---------- Gaussian ----
        Gauss[row, col] = exp(-s)
    else
        # ---------- SOAR ----
        SOAR[row, col] = 0
        Gauss[row, col] = 0
    end
    # if row >= 3
    #     break
    # end
    if row%1000 == 0 && c_lat == 1
        println("Row: $row / $n")
    end
end

##
println("Plotting...")
plotResults(SOAR, Gauss, nₕ, nᵥ)
println("Done!")

# Custom size correlation function with proper Rossbiy deformation taked into account.