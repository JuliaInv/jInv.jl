using jInv.Utils
using Test

v  = rand(1:200,54)
vu = unique(v)
w  = randn(100)

data = (vu,w)
# testing sortpermFast

for k=1:length(data)

	res1 = sortpermFast(data[k])
	res2 = sortperm(data[k])
	b2   = data[k][res2]

	@test all(res1[1] .== res2)
	@test all(res1[2] .== b2)
end

# Test two argument version
vin  = copy(v)
vuin = copy(vu)
d    = rand(1:1000,54)
du   = d[1:length(vu)]
din  = copy(d)
duin = copy(du)

ii,vs  = sortpermFast(v)
jj,vus = sortpermFast(vu)

vin,din   = sortpermFast(vin,din)
vuin,duin = sortpermFast(vuin,duin)

@test vin == vuin == vus
@test duin == du[jj]
@test length(din) == length(duin)
