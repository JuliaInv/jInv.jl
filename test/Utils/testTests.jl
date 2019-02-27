using jInv.Utils
using Test


# test checkDerivative
f = x-> x^2
dfRight = x -> 2*x
p,e,o = checkDerivative(f,dfRight,randn(),out=true)
@test p==true

F = x->x.^2
dFRight = x-> Diagonal(2*x)

x0 = randn(10).+2
p,e,o = checkDerivative(F,dFRight,x0,out=true)
@test p==true

dfWrong = x-> x
dFWrong = x->randn(length(x),length(x))
p,e,o = checkDerivative(f,dfWrong,randn(),out=false)
@test p==false

x0 = randn(10)
p,e,o = checkDerivative(F,dFWrong,x0,out=false)
@test p==false
