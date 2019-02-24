
using Test
using jInv.Mesh

h = rand(6)
msh = getTensorMesh3D(h,h,h)
m = rand(msh.nc)

Ne, Qe, activeEdges = getEdgeConstraints(msh);
Nf,Qf = getFaceConstraints(msh);

Curl = getCurlMatrix(msh)
Msig = getEdgeMassMatrix(msh, m)
Mmu = getFaceMassMatrix(msh, m)

CurlC = Qf  * Curl * Ne
MsigC = Ne' * Msig * Ne
MmuC  = Nf' * Mmu * Nf

@test Curl==CurlC
@test Msig==MsigC
@test Mmu==MmuC
