export removeZeros, complex2real, real2complex

#--------------------------------------------------------------------

function removeZeros( A::Array )
# Make array shorter by removing all zero values.
A = A[A .!= 0]
return A
end # function removeZeros

#---------------------------------------------------------------------

function complex2real( a::Array{Complex128} )
# Convert from [a1 a2 ...] to [ ar1 ai1 ar2 ai2 ... ]
n = length(a)
b = Array{Float64}(n*2)

j = 1
for i = 1:n
   b[j], b[j+1] = reim(a[i])
   j += 2
end  # i

return b   
end # function complex2real

#---------------------------------------------------------------------

function real2complex( b::Array{Float64} )
# Convert from [ br1 bi1 br2 bi2 ... ] to [b1 b2 ...]
n = length(b)
a = Array{Complex128}(div(n,2))

j = 1
for i = 1:2:n
   a[j] = complex( b[i], b[i+1] )
   j += 1
end  # i

return a
end  # function real2complex
