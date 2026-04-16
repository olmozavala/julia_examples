module MyModule
using Plots

mutable struct MyData
    x
    y
end

function simplePlot(data::MyData)
    print("The current backend in use is: ", backend())
    println(" $(data.x) and $(data.y)")
    ## Making some simple plots
    x = -2π:.1:3π
    y = sin.(x)

    plot(x,y)
    # plot([1 2 3 54  65])
end

export MyData, simplePlot
end