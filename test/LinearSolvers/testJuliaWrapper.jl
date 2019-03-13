using jInv.Mesh
using jInv.LinearSolvers
using SparseArrays
using LinearAlgebra
using Test
println("===  Test Julia wrapper with symmetric (real) matrix ====");
domain = [0.0, 1.0, 0.0, 1.0];
n      = [33,17];
Mr     = getRegularMesh(domain,n)
G      = getNodalGradientMatrix(Mr);
m      = sparse(Diagonal(exp.(randn(size(G,1)))));
Ar     = G'*m*G;
Ar     = Ar + 1e-1*norm(Ar,1)*sparse(1.0I,size(Ar,2),size(Ar,2));
N      = size(Ar,2);
b      = Ar*rand(N,1);
bs     = Ar*rand(N,3);
Bs     = (b,sparse(b),bs,sparse(bs));

sNonSym = getJuliaSolver(sym = 0);
sSym    = getJuliaSolver(sym = 1);
sGenSym = getJuliaSolver(sym = 2);
sols    = (sNonSym,sSym,sGenSym)

for j=1:length(Bs)
	print("nrhs=$(size(Bs[j],2)), issparse(rhs)=$(issparse(Bs[j])) : ")
	for k=1:length(sols)
		print("sym=$(sols[k].sym), ")
		xtt1,  = solveLinearSystem(Ar,Bs[j],sols[k]);
		@test norm(Ar*xtt1-Bs[j],Inf)/norm(Bs[j],Inf) < 1e-10
	end
	print("\n")
end
println("\n")


println("===  Test Julia wrapper with shifted (complex) Laplacian ====");
Ar     = G'*m*G +  (1+1im)*sparse(1.0I,size(Ar,2),size(Ar,2));
b      = Ar*rand(N,1);
bs     = Ar*rand(N,3);
Bs     = (b,sparse(b),bs,sparse(bs));

sNonSym = getJuliaSolver(sym = 0);
sSym    = getJuliaSolver(sym = 1);
sGenSym = getJuliaSolver(sym = 2);
sols    = (sNonSym,sSym,sGenSym)

for j=1:length(Bs)
	print("nrhs=$(size(Bs[j],2)), issparse(rhs)=$(issparse(Bs[j])) : ")
	for k=1:length(sols)
		print("sym=$(sols[k].sym), ")
		xtt2,  = solveLinearSystem(Ar,Bs[j],sols[k]);
		@test norm(Ar*xtt2-Bs[j],Inf)/norm(Bs[j],Inf) < 1e-10
	end
	print("\n")
end
println("\n")


println("===  Test Julia Wrapper: nonsymmetric matrices ====");
n = 50
sNonSym  = getJuliaSolver(sym = 0);
sNonSym  = copySolver(sNonSym);
A = sprandn(n,n,5/n) + 10*sparse(1.0I,n,n)
B = randn(n)

X, = solveLinearSystem(A,B,sNonSym,0);
@test norm(A*X - B)/norm(B) < 1e-10

# multiple sparse rhs - 
# Bs = sprandn(n,10,.3)
# X, = solveLinearSystem(A,Bs,sNonSym,0);
# @test norm(A*X - Bs,Inf)/norm(Bs,Inf) < 1e-10


X, = solveLinearSystem(A,B,sNonSym,1);
@test norm(A'*X - B)/norm(B) < 1e-10

X, = solveLinearSystem(A,B,sNonSym,1);
@test norm(A'*X - B)/norm(B) < 1e-10

X, = solveLinearSystem(A,B,sNonSym,0);
@test norm(A*X - B)/norm(B) < 1e-10

println("===  End Test Julia Wrapper ====");
