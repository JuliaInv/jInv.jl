using jInv.Mesh
using Test

# setup 3D problem
nc = [5, 7, 8]
x0 = [0.1,0.2, 0.4]
domain = [x0[1], 4, x0[2], 2, x0[3], 6]
h   = (domain[2:2:end]-domain[1:2:end])./nc
h1  = h[1]*ones(nc[1])
h2  = h[2]*ones(nc[2])
h3  = h[3]*ones(nc[3])

Mt = getTensorMesh3D(h1,h2,h3,x0)
Mr = getRegularMesh(domain,nc)

print("\ttest cell-centered axes...")
xt = getCellCenteredAxes(Mt)
xr = getCellCenteredAxes(Mr)
for d=1:3
	@test norm(xt[d]-xr[d],Inf)<1e-14
end
print("passed\n")

print("\ttest nodal axes...")
xt = getNodalAxes(Mt)
xr = getNodalAxes(Mr)
for d=1:3
	@test norm(xt[d]-xr[d],Inf)<1e-14
end
print("passed\n")

print("\ttest cell-centered grid...")
xt = getCellCenteredGrid(Mt)
xr = getCellCenteredGrid(Mr)
@test norm(xt-xr,Inf)<1e-14
print("passed\n")

print("\ttest nodal grid...")
xt = getNodalGrid(Mt)
xr = getNodalGrid(Mr)
@test norm(xt-xr,Inf)<1e-14
print("passed\n")

print("\ttest face grid...")
xt = getFaceGrids(Mt)
xr = getFaceGrids(Mr)
for d=1:3
	@test norm(xt[d]-xr[d],Inf)<1e-14
end
print("passed\n")

print("\ttest edge grid...")
xt = getEdgeGrids(Mt)
xr = getEdgeGrids(Mr)
for d=1:3
	@test norm(xt[d]-xr[d],Inf)<1e-14
end
print("passed\n")

print("\ttest volume...")
xt = getVolume(Mt)
xr = getVolume(Mr)
@test norm(xt-xr,Inf)<1e-14
xt = getVolumeInv(Mt)
xr = getVolumeInv(Mr)
@test norm(xt-xr,Inf)<1e-14
print("passed\n")

print("\ttest area...")
xt = getFaceArea(Mt)
xr = getFaceArea(Mr)
@test norm(xt-xr,Inf)<1e-14
xt = getFaceAreaInv(Mt)
xr = getFaceAreaInv(Mr)
@test norm(xt-xr,Inf)<1e-14
print("passed\n")

print("\ttest length...")
# test length
xt = getLength(Mt)
xr = getLength(Mr)
@test norm(xt-xr,Inf)<1e-14
# test length inverse
xt = getLengthInv(Mt)
xr = getLengthInv(Mr)
@test norm(xt-xr,Inf)<1e-14
print("passed\n")

print("\ttest counting...")
@test Mt.nc==Mr.nc
@test Mt.ne==Mr.ne
@test Mt.nf==Mr.nf
@test Mt.n ==Mr.n
print("passed\n")

print("\ttest nodal gradient matrix...")
Gr  = getNodalGradientMatrix(Mr)
Gt  = getNodalGradientMatrix(Mt)
@test norm(Gr-Gt,1)/norm(Gr,1) < 1e-12
print("passed\n")

print("\ttest divergence matrix...")
Dr  = getDivergenceMatrix(Mr)
Dt  = getDivergenceMatrix(Mt)
@test norm(Dr-Dt,1)/norm(Dr,1) < 1e-12
print("passed\n")

print("\ttest curl matrix...")
Cr  = getCurlMatrix(Mr)
Ct  = getCurlMatrix(Mt)
@test norm(Cr-Ct,1)/norm(Cr,1) < 1e-12
print("passed\n")

print("\ttest face average matrix...")
Afr = getFaceAverageMatrix(Mr)
Aft = getFaceAverageMatrix(Mt)
@test norm(Afr-Aft,1)/norm(Afr,1) < 1e-12
print("passed\n")

print("\ttest edge average matrix...")
Aer = getEdgeAverageMatrix(Mr)
Aet = getEdgeAverageMatrix(Mt)
@test norm(Aer-Aet,1)/norm(Aer,1) < 1e-12
print("passed\n")
