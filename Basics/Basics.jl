using LinearAlgebra
using Printf
using Plots
##

α = 3 # Unicode varibles
:a  # Defines a symbol that you can use later 
println(typeof(:a))

## --------- Exceptions --------
x = 1
println("--------- Exceptions --------");
try 
   sqrt("sopas")
catch e
   println("Error: $e")
finally
   x = 3
end
println("x = $x")

## ------- Struct ------
mutable struct Point
   x
   y
end
x = Point(3,4)
print(x)

##--------------- Maps, filters, lists comprehension (Apply a function to a vector)
println(" #--------------- Maps (Apply a function to a vector)")
arr = [1,2,3,5]
byThen = map(x->x*10, arr)
println(byThen);# Uses the multiply by 10 map
println(filter(x -> iseven(x),arr));# Filter the array  by only evens
println(Float64[x^2 for x in 1:10])# Initializes an array with a function (first 10 squares)
println([x^2 for x in 1:10])# Initializes an array with a function (first 10 squares (ints))
# a = [true true true]
a = [true true true]
@show reduce(&, a)

## ---------- Matrices --------
println("# ---------- Matrices --------")
mat = [1 0; 0 5]
println(mat)
b = [10 432]
#println(mat/b) # Solving Ax=b

## ---------- Dictionaries--------
println("# ---------- Dictionaries--------")
# Basic dict
d = Dict("sopas"=>31, "perico"=>123)  # Same type
println("Sopas: ", d["sopas"])
for k in keys(d)
   println(k)
end
d = Dict(1=>rand(3,2))  # Same type
println(typeof(d))
println(keys(d))
println(d[1])

## Similar with a named tuple
d = (sopas = 31, per = 5)
d[:sopas]
d[:per]

## With symbols
dog = Dict(:name => "hiro", :age => "old")
println(dog[:name])

## From a different file 
println("Loading external modules")
include("MyModule.jl")
data = MyModule.MyData([1 2 3], [4 5 6])
println(data)
println(data.x)
MyModule.simplePlot(data)

## ------- Macros --------
macro sayhi(name)
   println("Yeah babe!!!!! You are $name")
end
@sayhi "Olmuz"

