
using jInv.InverseSolve
using jInv.Utils

n = 200
x = rand(n) + im*rand(n)
A = sprand(n,2*n,0.1)

Hx = A[1:2:end,1:2:end]*real(x) +
     A[1:2:end,2:2:end]*imag(x) +
     im*A[2:2:end,1:2:end]*real(x) +
     im*A[2:2:end,2:2:end]*imag(x)

@test all(abs.(HessMatVec(A,x)-Hx).<1e-12)
