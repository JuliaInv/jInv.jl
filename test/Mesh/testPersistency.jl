using jInv.Mesh
using jInv.Utils
using Test

# setup 3D mesh
nc = [5, 7, 8]
x0 = rand(3)
domain = [x0[1], 4, x0[2], 2, x0[3], 6]
h   = (domain[2:2:end]-domain[1:2:end])./nc
h1  = h[1]*ones(nc[1])
h2  = h[2]*ones(nc[2])
h3  = h[3]*ones(nc[3])

Mt = getTensorMesh3D(h1,h2,h3,x0)
Mt2 = getTensorMesh3D(h1,h2,h3,x0)
Mr = getRegularMesh(domain,nc)
Mr2 = getRegularMesh(domain,nc)

# test if Mt and Mt2 are equal
@test Mt == Mt2
@test Mr == Mr2

Meshes = (Mt,Mt2,Mr,Mr2)

function buildOps!(M)
	D = getDivergenceMatrix(M)
	G = getNodalGradientMatrix(M)
	C = getCurlMatrix(M)
	Af = getFaceAverageMatrix(M)
	Ae = getEdgeAverageMatrix(M)
	An = getNodalAverageMatrix(M)
	V  = getVolume(M)
	F  = getFaceArea(M)
	L  = getLength(M)
	Vi = getVolumeInv(M)
	Fi = getFaceAreaInv(M)
	Li = getLengthInv(M)
	Lap = getNodalLaplacianMatrix(M)
end

# first all fields should be empty
for k=length(Meshes)
	@test isempty(Meshes[k].Div )
	@test isempty(Meshes[k].Grad)
	@test isempty(Meshes[k].Curl)
	@test isempty(Meshes[k].Af  )
	@test isempty(Meshes[k].Ae  )
	@test isempty(Meshes[k].An  )
	@test isempty(Meshes[k].V   )
	@test isempty(Meshes[k].F   )
	@test isempty(Meshes[k].L   )
	@test isempty(Meshes[k].Vi  )
	@test isempty(Meshes[k].Fi  )
	@test isempty(Meshes[k].Li  )
	@test isempty(Meshes[k].nLap)
end

# build all matrices
buildOps!(Mt)
buildOps!(Mt2)
buildOps!(Mr)
buildOps!(Mr2)
@test Mt == Mt2
@test Mr == Mr2


# now all fields should be non-empty
for k=length(Meshes)
	@test !isempty(Meshes[k].Div )
	@test !isempty(Meshes[k].Grad)
	@test !isempty(Meshes[k].Curl)
	@test !isempty(Meshes[k].Af  )
	@test !isempty(Meshes[k].Ae  )
	@test !isempty(Meshes[k].An  )
	@test !isempty(Meshes[k].V   )
	@test !isempty(Meshes[k].F   )
	@test !isempty(Meshes[k].L   )
	@test !isempty(Meshes[k].Vi  )
	@test !isempty(Meshes[k].Fi  )
	@test !isempty(Meshes[k].Li  )
	@test !isempty(Meshes[k].nLap)
end

clear!(Mt)
clear!(Mt2)
clear!(Mr)
clear!(Mr2)

@test Mt == Mt2
@test Mr == Mr2

# now all matrices should be gone
for k=length(Meshes)
	@test isempty(Meshes[k].Div )
	@test isempty(Meshes[k].Grad)
	@test isempty(Meshes[k].Curl)
	@test isempty(Meshes[k].Af  )
	@test isempty(Meshes[k].Ae  )
	@test isempty(Meshes[k].An  )
	@test isempty(Meshes[k].V   )
	@test isempty(Meshes[k].F   )
	@test isempty(Meshes[k].L   )
	@test isempty(Meshes[k].Vi  )
	@test isempty(Meshes[k].Fi  )
	@test isempty(Meshes[k].Li  )
	@test isempty(Meshes[k].nLap)
end
