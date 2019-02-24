# Testing sparseUtils
using jInv.Utils
using Test

# test sdiag
a = randn(13)
@test a == diag(sdiag(a))
