using jInv.Utils
using Test


a = rand(1:300,300)
ua = unique(a)

# testing uniqueidx
b,ii,jj = uniqueidx(a)
b2      = sortunique(a)
t1 = sort(ua)
t2 = sortperm(ua)

@test all(b .== t1)
@test all(a[ii] .==ua[t2])
@test all(b[jj].==a)
@test b == b2
