
using jInv.Mesh

h = rand(6)
msh = getTensorMesh3D(h,h,h)

Ne, Qe, activeEdges = getEdgeConstraints(msh);
Nf,Qf = getFaceConstraints(msh);

Curl = getCurlMatrix(msh)
Curl = Qf  * Curl * Ne
Msig = Ne' * Diagonal(ones(msh.ne)) * Ne
Mmu  = Nf' * Diagonal(ones(msh.nf))  * Nf
