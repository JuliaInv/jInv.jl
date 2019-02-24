using jInv.Mesh
using Test

# setup 3D problem
nc = [5, 7, 8]
x0 = rand(3)
domain = [x0[1], 4, x0[2], 2, x0[3], 6]
h   = (domain[2:2:end]-domain[1:2:end])./nc
h1  = h[1]*ones(nc[1])
h2  = h[2]*ones(nc[2])
h3  = h[3]*ones(nc[3])

Mt = getTensorMesh3D(h1,h2,h3,x0)
Mr = getRegularMesh(domain,nc)
Mt2 = getTensorMesh3D(h1+rand(nc[1]),h2+rand(nc[2]),h3+rand(nc[3]),x0)

Meshes = (Mt,Mr,Mt2)

for k=1:length(Meshes)
	M = Meshes[k]
	print("\ttesting differential operators for $(typeof(M))...")
	@test norm(getDivergenceMatrix(M)*getCurlMatrix(M),Inf)<1e-12
	
	Lap = getNodalLaplacianMatrix(M)
	G   = getNodalGradientMatrix(M)
	@test norm(Lap*ones(size(Lap,2))/M.nc) < 1e-10
	@test norm(G'*G- Lap,Inf) <1e-12	
	print("passed!\n")
end


