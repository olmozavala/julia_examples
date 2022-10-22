using Printf

##  ------ Functions
function f(x) # General function
   "1_$(x)"
end
y = ["a" "b" "c"]
println("Element wise function call: ", f.(y)) # Calling it and Apply element wise

## ------- Parameters options (type, default, named, additional)
# Named variables after a semicolon
function f(a::Int, b::Int=3, x...; name="Sopas") 
   println("a=$a b=$b name=$name x=$x")
   println("Len x = $(length(x))")
end
println("Calling function with: f( name='Yeah', 1, 34, 55,2,2)")
# f(1.5,34,55,2,2) # This will generate a type error
f(name="Yeah", 1,34, 55,2,2)

println("\nCalling function with: f(1)")
f(1)

## ------ Array in function params
function farr(x::Array) # General function
end

## ------ Return Parameters
println(" # ------ Return ------")
function myfunc(x,y)
   return x + y, x*y; #Return more than one value
end
out = myfunc(23,41)
println("outputs: ", out)

## ------ Inline and anonymous functions
inlineFunc(x,y) = x*y
println("Inline function: ", inlineFunc(3,2)); 
x = (n) -> n*n; # Example of anonymous function
println("Anonymous function:", x(3)); # Concatenates the array with the other params

## ------- Functions of functions
function f_of_f(f,x,y)
   return f(x,y)
end
println("Function of function: f_of_f(inlineFunc,3,3)) = ", f_of_f(inlineFunc,3,3)); # Concatenates the array with the other params

## ------- Macros --------
macro sayhi(name)
   println("Yeah babe!!!!! You are $name")
end
@sayhi "Olmuz"

## ------- External Functions --------
include("MyFunctionsFile.jl") 
println(addPlusTree(2,2))