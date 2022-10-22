##  ------ if
if false
   print("true")
elseif true
   print("false")
else
   print("nop")
end
##  ------ Inline if
println(3 < 5 ? "Yes" : "no")

## -------- For
println("For: ")
for i = 1:5
   println(i)
end

# for i in ["sopas", "perico", "lico"]
for i ∈ ["sopas", "perico", "lico"]
   println(i)
end

## ennumerate
for (i,val) ∈ enumerate(["sopas", "perico", "lico"])
   println("$i and value: $val")
end

## Nested loop
println("Nested loop:")
for i = 1:2, j = 3:4
   println((i, j))
end

## ----- While
print("While")
i=1
while i <= 5
   println(i)
   i += 1
end

##  ------ Comprehensions
println("Comprehensions")
A = [x^2 for x in 1:3]
println([x^2  for x in 1:4 if x > 2])
println(A)

