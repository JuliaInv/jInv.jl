using Test

A    = sprandn(100,10,.1)
pFor = LSparam(A,[])
m0   = randn(10)

# checkDerivative for ForwardProbTypes
passedDerivative, = checkDerivative(m0,pFor)
@test passedDerivative

# adjoint test
passedAdjoint, = adjointTest(m0,pFor)
@test passedAdjoint
