# Helper function to make orbit plot
function plot_od(odiag, pvalues, n = 1000)
    println("Ploting...")
    L = length(pvalues)
    x = Matrix{Float64}(undef, n, L)
    y = Matrix{Float64}(undef, n, L)
    for j in 1:L
        x[:,j] .= pvalues[j]
        y[:,j] .= odiag[j]
    end
    scatter(x, y, 
        markersize=0.1, markeralpha = 0.3, markercolor="black",
        leg=false, title="Bifurcation graph", 
        html_output_format=:png, size=(1000,500))
end