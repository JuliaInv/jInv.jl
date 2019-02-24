using jInv.Utils
using Test

# testing variousUtils

A = rand(0.:10., 10) + im*rand(0.:10., 10)

@test real(A) == complex2real(A)[1:2:end]
@test imag(A) == complex2real(A)[2:2:end]
@test A == real2complex(complex2real(A))
