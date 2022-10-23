# https://nextjournal.com/sosiris-de/pde-2018

using Plots
gr()

# ===================== Solving linear Poisson (1D) uₓₓ = sin(2πx) =====================
Δx = 0.01
x = Δx:Δx:1-Δx # Solve only for the interior: the endpoints are known to be zero!
N = length(x)
B = sin.(2π*x)

# We solve uₓₓ with central difference u_xx = u_{i+1} - 2u_i + u_{x-1} / Δx²
# If u(0) = u(1) = 0 then 
A = zeros(N,N)
for i in 1:N, j in 1:N
  abs(i-j)<=1 && (A[i,j]+=1)
  i==j && (A[i,j]-=3)
end
A = A/(Δx^2)

# Now we want to solve AU=B, so we use backslash:
U = A\B
plot(x, -sin.(2π*x)./(4π^2), label="True")
scatter!([0;x;1], [0;U;0])

## ===================== Solving linear Poisson (2D) uₓₓ + u_yy = sin(2πx) + sin(2πy) =====================
Δx = 0.01
Δy = 0.01
x = Δx:Δx:1-Δx # Solve only for the interior: the endpoints are known to be zero!
y = Δy:Δy:1-Δy # Solve only for the interior: the endpoints are known to be zero!
N = length(x)
M = length(y)

B = sin.(2π*x)' .* ones(size(x)) .+ 
    transpose(sin.(2π*y)' .* ones(size(y)))


# We solve uₓₓ with central difference u_xx = u_{i+1} - 2u_i + u_{x-1} / Δx²
# If u(0) = u(1) = 0 then 
A = zeros(N,N)
for i in 1:N, j in 1:N
  abs(i-j)<=1 && (A[i,j]+=1)
  i==j && (A[i,j]-=3)
end
A = A/(Δx^2)

# Now we want to solve AU=B, so we use backslash:
U = A\B
plot(x, -sin.(2π*x)./(4π^2), label="True")
scatter!([0;x;1], [0;U;0])


# ===================== Solving nonlinear Poisson uₓₓ = f(u) =====================