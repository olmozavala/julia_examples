using LinearAlgebra
using Printf

## ------- Arrays ------
# https://docs.julialang.org/en/v1/manual/arrays/#
println("------- Arrays ------")
arr = Float64[]
push!(arr, 123.2) # Push into the array
println("Array: ", arr);
rx = rand(10)
println("Random array: ", rx)
println("Random array: $([@sprintf("%.2f",x) for x in rx])")

##
println("------- Properties------")
arr_col1 = ones(3); # Col vector
arr_col2 = [2,3,5] # Col vector
arr_row = [1 1 1] # Row vector
println("1 arr_col1: ", arr_col1)
println("Reshape: ", reshape(arr_col1, (3,1)))
println("2 arr_col2: ", arr_col2)
println("3 arr_row: ", arr_row)
println("1 Size: ", size(arr_col1))
println("2 Size: ", size(arr_col2))
println("3 Size: ", size(arr_row))
println("Length: ", length(arr_row))
println("Concatenation: ", [arr_col1..., arr_col2...])
println("Array mult by element: ", arr.*arr_col2)
println("Array mult (dot product): ",dot(arr,arr_col2))
println("Array mult matrix (row by col): ",transpose(arr_col1)*arr_col2)
println("Array mult matrix (row by col): ",arr_row*arr_col2)
println("Array mult matrix (col by row): ",arr_col2*arr_row)

##  ------ Basic methods
println(arr_col1.min)

##  ------ Index
x = collect(0:1:10) # An equivalent to range
println("x ", x)
println("x[1:5] → ", x[1:5])
println("x[1:2:5] →  ", x[1:2:5])
println("x[end] →  ", x[end])
println("x[1:2:end] →  ", x[1:2:end])
println("x[1:2:end-5] →  ", x[1:2:end-5])
println("circshift(x,2) → ", circshift(x,2))

## ---------- Matrices --------
println("# ---------- Matrices --------")
println(ones(5))
println(zeros(3,2))
println(size(rand(3,5)))
mat = [1 0; 0 5]
println(mat)
b = [10 432]
println(mat/b) # Solving Ax=b

## ---------- Tuples (inmmutable)
tup = (1, "sopas")
println(tup)
# tup[1] = 3 # WIll fail