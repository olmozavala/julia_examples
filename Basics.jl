using LinearAlgebra
using Printf
##

# Variable declaration
x = 9
str = "This is a string, var value: $x"
println(str)

α = 3 # Unicode varibles
:a  # Defines a symbol that you can use later 
typeof(:a)

#Evaluating variables and expressions in a string  Strings
println("Evaluating variables and expressions in a string")
println("Testing value of x = $x")
println("Testing value of x = $(1.423)")
println("Testing value of x+2 = $(x+2)")
println("Leading zeros = $(lpad(5,3,'0'))")
@printf("Number of decimals= %.2f",rand(1)[1])

## ---- Testing regex -------
println(" ---- Testing regex -------");
myregex = r".+@.+"
println("Does it match: ", occursin(myregex,"so@pas"))
println("Does it match: ", occursin(myregex,"sopas"))


## ------- Struct ------
mutable struct Point
   x
   y
end
x = Point(3,4)
print(x)

## ------- Arrays ------
# https://docs.julialang.org/en/v1/manual/arrays/#
println("------- Arrays ------")
arr = Float64[]
push!(arr, 123.2)
println("Array: ", arr);
println( rand(10))

arr = ones(3);
secarr = [2,3,5];
println("Size: ", size(secarr))
println("Length: ", length(secarr))
println("Array mult: ",arr.*secarr)
println("Array mult: ",dot(arr,secarr))
println("Array mult: ",transpose(arr)*secarr)
##  ------ Index
x = collect(0:1:10) # An equivalent to range
println(x)
println(x[1:5])
println(x[1:2:5])
println(x[end])
println(x[1:2:end])
println(circshift(x,2))

##  ------ Comprehensions
println("Comprehensions")
A = [x^2 for x in 1:3]
display([x^2  for x in 1:3 if x > 2])
println(A)

##  ------ Functions
function f(x) # Basic
   string(x , "1")
end
y = ["a" "b" "c"]
println(f.(y))
function f(a,b=3,x...) # Varagrs
   println(b)
   println(x)
   println(length(x))
end
display(f(1,34,55,2,2))


include("MyFunctionsFile.jl")
println(addPlusTree(2,2))

##  ------ if
if false
   print("true")
else
   print("false")
end

## Control flow
print("Loops")
i =1
while i <= 5
   println(i)
   global i += 1
end

for i = 1:5
   println(i)
end

# for i in ["sopas", "perico", "lico"]
for i ∈ ["sopas", "perico", "lico"]
   println(i)
end

# nested loop
for i = 1:2, j = 3:4
   println((i, j))
end


## ------ Functions ------
println(" # ------ Functions ------")
function myfunc(x,y)
   return x + y, x*y; #Return more than one value
end
println("Simple function with two outputs: ", myfunc(23,41));

function myfuncVariable(args...)
   for i in args
       print(i,",")
   end
   return 1
end
println("Multiple arguments: ", myfuncVariable(1,2,[3,4]))
println("Multiple arguments: ", myfuncVariable(1,2,[3,4]...)); # Concatenates the array with the other params

function defVal(x, y=3)
   return x+y
end
println("Default value: ", defVal(1)); # Concatenates the array with the other params
println("Default value: ", defVal(1, 8)); # Concatenates the array with the other params

inlineFunc(x,y) = x*y
println("Inline function: ", inlineFunc(1234,123)); # Concatenates the array with the other params

function funcOffunc(f,x,y)
   return f(x,y)
end
println("Function of function: ", funcOffunc(inlineFunc,3,3)); # Concatenates the array with the other params

x = (n) -> n*n; # Example of anonymous function
println("Anonymous function:", x(3)); # Concatenates the array with the other params

# From a different file 
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


## --------- Exceptions --------
println("--------- Exceptions --------");
try 
   sqrt("sopas")
catch e
   println("Error: $e")
end

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

# Similar with a named tuple
d = (sopas = 31, per = 5)
d[:sopas]
d[:per]

# With symbols
dog = Dict(:name => "hiro", :age => "old")
println(dog[:name])

