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
M2D = getRegularMesh(domain[1:4],nc[1:2])
Mt2 = getTensorMesh3D(h1+rand(nc[1]),h2+rand(nc[2]),h3+rand(nc[3]),x0)

Meshes = (Mt,Mr,Mt2,M2D)

for k=1:length(Meshes)
	M = Meshes[k]
	print("\ttesting nodal averaging for $(typeof(M))...")
	Av = getNodalAverageMatrix(M)
	xn = getNodalGrid(M)
	xc = getCellCenteredGrid(M)
	@test norm(Av*xn-xc,Inf)<1e-12
	
	print("passed!\n")
end



