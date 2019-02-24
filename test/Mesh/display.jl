using jInv.Mesh
using Test

# regular mesh
domain = [1.2 2.4 0.0 1.1 0.0 2.3]
n      = [3 4 5]
M3D      = getRegularMesh(domain,n)
display(M3D)
M2D      = getRegularMesh(domain[1:4],n[1:2])
display(M2D)

# tensor mesh
h1 = rand(4)
h2 = rand(5)
h3 = rand(6)
MT = getTensorMesh3D(h1,h2,h3)
display(MT)