using Plots
gr()

function getGaussShape(n::Int, σ::Float64, d::Int = 1, μ::Float64=0.0)
    x = getGauss(n, σ, d, μ)
    x = x./(maximum(x))
    return x
end

function getGauss(n::Int, σ::Float64, d::Int = 1, μ::Float64=0.0)
    # TODO this is not true in 2D
    # Simple function to obtain a Gaussian function
    # g(x) = (1/(σ*√(2π))).*exp.(-.5((x.-μ).^2)/σ^2)
    g(x) = (1/(σ*√(2π))).*exp.(-.5((x.-μ).^2)/σ^2)
    if d == 1
        # In this case it is 1D
        println(1)
        x = range(-n/2, n/2, length=n)
        return g(x) 
    end
    if d == 2
        if μ == 0.0
            μ = (0.0, 0.0)
        end
        A = zeros(n,n)
        # In this case it is 2D
        x = range(-n/2, n/2, length=n)
        # A = (1/(σ*√(2π))).*exp.(-.5(x.-μ)/σ^2)
        for (i, y) in enumerate(x)
            A[i,:] = exp.(-.5((y-μ[1])^2 .+ (x.-μ[2]).^2)/σ^4)
        end
        return A
    end
end

# scatter(getGauss(100, 10.0))
# heatmap(getGauss(100, 3.0, 2))
# heatmap(getGauss(100, 3.0, 2, (0.0,10.5)))