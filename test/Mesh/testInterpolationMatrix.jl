using LinearAlgebra
using jInv.Mesh
using Test


print("   test getInterpolationMatrix (RegularMesh) ... ")
domain = zeros(6); domain[2:2:end] = rand(3).+1; domain[1:2:end] = -rand(3).-1
n1    = [2; 3; 4]
n2    = 2*n1
for dim=2:3
	print(" dim=",dim)
	Mc    = getRegularMesh(domain[1:2*dim],n1[1:dim])
	Mf    = getRegularMesh(domain[1:2*dim],n2[1:dim])
	Pcf   = getInterpolationMatrix(Mc,Mf)
	Pfc   = getInterpolationMatrix(Mf,Mc)
	Xc    = getCellCenteredGrid(Mc)
	Xf    = getCellCenteredGrid(Mf)

	@test isapprox(sum(Pfc,dims=2),ones(size(Pfc,1)))
	@test norm(Xc-Pfc*Xf)/norm(Xc) < 1e-13
end
print(" passed\n")
print("   test getInterpolationMatrix (TensorMesh) ... ")
x0    = randn(3)
h1    = rand(4)
h2    = rand(6)
h3    = rand(8)
H1    = h1[1:2:end]+h1[2:2:end]
H2    = h2[1:2:end]+h2[2:2:end]
H3    = h3[1:2:end]+h3[2:2:end]

Mct    = getTensorMesh3D(H1,H2,H3,x0)
Mft    = getTensorMesh3D(h1,h2,h3,x0)
Pcft   = getInterpolationMatrix(Mct,Mft)
Pfct   = getInterpolationMatrix(Mft,Mct)
Xc    = getCellCenteredGrid(Mct)
Xf    = getCellCenteredGrid(Mft)

@test isapprox(sum(Pfct,dims=2),ones(size(Pfct,1)))
@test norm(Xc-Pfct*Xf)/norm(Xc) < 1e-13
print("passed\n")
