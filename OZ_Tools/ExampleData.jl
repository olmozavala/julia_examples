using Plots
gr()

function getGauss(n::Int, σ::Float64, d::Int = 1)
    # TODO this is not true in 2D
    # Simple function to obtain a Gaussian function
    μ = 0
    g(x) = (1/(σ*√(2π))).*exp.(-.5((x.-μ).^2)/σ^2)
    if d == 1
        # In this case it is 1D
        println(1)
        x = range(-n/2, n/2, length=n)
        return g(x) 
    end
    if d == 2
        A = zeros(n,n)
        # In this case it is 2D
        x = range(-n/2, n/2, length=n)
        # A = (1/(σ*√(2π))).*exp.(-.5(x.-μ)/σ^2)
        for (i, y) in enumerate(x)
            A[i,:] = g(y)*(1/(σ*√(2π))).*exp.(-.5((x.-μ).^2)/σ^2)
        end
        return A
    end
end

scatter(getGauss(100, 10))
# heatmap(getGauss(100, 5, 2))