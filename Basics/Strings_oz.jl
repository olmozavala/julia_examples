using Printf
##

x = 1.2425324253
println("------------- Strings ---------------------")
println("Testing value of x = $x")
println("Testing value of $(1.423)")
println("Leading zeros = $(lpad(5,3,'0'))")
println("""
        Here some unstructured Strings
        Yeah babe! """)

## ---- Printf
# % c,s,d,f, .nf, e
@printf("Number of decimals= %.2f \n", x)
@printf("Testing value of x+2 = %.2f \n", (x+2))

## ---- Testing regex -------
println(" ---- Testing regex -------")
myregex = r".+@.+"
println("Does it match: ", occursin(myregex,"so@pas"))
println("Does it match: ", occursin(myregex,"sopas"))


## ---- Styled
printstyled("Yeah babe", bold=true, color=:blue)