using   ModelingToolkit
# https://mtk.sciml.ai/stable/tutorials/symbolic_functions/

## Hello World (variables)
@variables x y  # Create some symbolic variables

z = x^2 + y

println("z: ", z)
println("Continue:", (z^2) * x + y^2)

## Array of variables
@variables u[1:3]  # Create an array of symbolic variables
println(u)

# ======== Derivatives ==========
Dx = Differential(x)
Dy = Differential(y)
println("Derivative of z:", Dx(z))
println("zₓ expanded: ", expand_derivatives(Dx(z)))#
println("z_y expanded: ", expand_derivatives(Dy(z)))#