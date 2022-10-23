using Plots
using Random
using ForwardDiff
# plotly()
gr()

## Create a syntetic 2D vector field. 
# A gyre centered at 0,0

dh = 1 
d =  π
x = -d:dh:d

X = (x * ones(length(x))')'
Y = x * ones(length(x))'
svec = .2
# Single gyre
# U = Y*svec
# V = -X*svec
# Source 
U = X*svec
V = Y*svec
# Sinc
# U = -X*svec
# V = -Y*svec

l = @layout[a b c d]
p1 = heatmap(X, title="X")
p2 = heatmap(Y, title="Y")
p3 = heatmap(U, title="U")
p4 = heatmap(V, title="V")
plot(p1, p2, p3, p4, layout=l, colorbar=:none)
## 
function plotvector(x,y,u,v,subsample)
    l = 1:length(vec(X))
    lt = shuffle(l)
    ls = lt[1:subsample:end]
    quiver(vec(X)[ls],vec(Y)[ls],quiver=(vec(U)[ls],vec(V)[ls]))
end

plotvector(X,Y,U,V,1)

## Divergence
d = (V[2:end,1:end-1] - V[1:end-1,1:end-1]) + (U[1:end-1,2:end] - U[1:end-1,1:end-1])

ld = @layout[a b]
# p1 = plot(x, U[1,:])
# =2 = plot( (U[1,2:end] - U[1,1:end-1])./dh)
p1 = plotvector(X,Y,U,V,1)
function div(U,V)
    m = size(U)[1]
    n = size(U)[2]
    d = zeros(m-1,n-1)
    for i=2:m-1
        for j=2:n-1
            d[i-1,j-1] = U[i,j+1] - U[i,j-1] + V[i+1,j] - V[i-1,j]
        end
    end
    return d
end
p2 = heatmap(div(V,U))
# a = V[1:4, 1:4]
# b = V[1:4, 1:4]
plot(p1, p2, layoUt=ld)


# DiffUsion
# Advection