# Testing sparseUtils
using jInv.Utils
using Base.Test

# test sdiag
a = randn(13)
@test a == diag(sdiag(a))
