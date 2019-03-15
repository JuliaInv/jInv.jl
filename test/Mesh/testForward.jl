
using jInv.Mesh
using Test

# setup 2D problem for Nodal DivSigGrad
nc = [5, 7]
x0 = rand(2)
domain = [x0[1], 4, x0[2], 2]
M = getRegularMesh(domain,nc)
print("\ttesting forward operators for $(typeof(M))...")
Lap = getNodalLaplacianMatrix(M)
Lsig = getNodalDivSigGradMatrix(M,ones(prod(M.nc)); saveMat = false, avN2C = avN2C_Nearest);
@test norm(Lsig- Lap,Inf) <1e-12	
print("passed!\n")


# setup 3D problem for Nodal DivSigGrad
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
	Mk = Meshes[k]
	print("\ttesting forward operators for $(typeof(M))...")
	Lap_k = getNodalLaplacianMatrix(M)
	Lsig_k = getNodalDivSigGradMatrix(M,ones(prod(M.nc)); saveMat = false, avN2C = avN2C_Nearest);
	@test norm(Lsig_k- Lap_k,Inf) <1e-12	
	print("passed!\n")
end


# setup 2D problem for Face Elasticity

print("\ttesting 2D elasticity operators for $(typeof(M))...")
nc = [5, 7]
x0 = rand(2)
domain = [x0[1], 4, x0[2], 2]
M = getRegularMesh(domain,nc)
mu = ones(prod(M.n));
lambda = ones(prod(M.n));
L = GetLinearElasticityOperator(M, mu,lambda);
L2 = GetLinearElasticityOperatorFullStrain(M, mu,lambda);
#println("Diff between two elasticity ops: ",norm(L-L2,1));
@test norm(L- L2,Inf) <1e-8	
print("passed!\n")


# setup 3D problem for Face Elasticity
print("\ttesting 3D elasticity operators for $(typeof(M))...")
n = [5,7,8];
domain = [0.0,1.0,0.0,2.0,0.0,1.0];
M = getRegularMesh(domain,n);
mu = ones(prod(M.n));
lambda = ones(prod(M.n));
L = GetLinearElasticityOperator(M, mu,lambda)
L2 = GetLinearElasticityOperatorFullStrain(M, mu,lambda)

# println("Diff between two elasticity ops: ",norm(L-L2,1));
@test norm(L- L2,Inf) <1e-8	
print("passed!\n")