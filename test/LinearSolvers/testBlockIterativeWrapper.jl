using jInv.Mesh
using jInv.LinearSolvers
using Test
using KrylovMethods

println("===  Example 2D DivSigGrad ====");

domain = [0.0, 1.0, 0.0, 1.0];
n      = [50,50];
Mr     = getRegularMesh(domain,n)
G      = getNodalGradientMatrix(Mr);
m      = spdiagm(exp.(randn(size(G,1))));
Ar     = G'*m*G;
Ar     = Ar + 1e-1*norm(Ar,1)*speye(size(Ar,2));
N      = size(Ar,2);
B      = Ar*rand(N,4);


sPCG   = getBlockIterativeSolver(KrylovMethods.blockCG,out=1,sym=1);
X,  = solveLinearSystem(Ar,B,sPCG);
@test vecnorm(Ar*X-B)/vecnorm(B) < sPCG.tol


IterMethod = (A,B;M=M,X=X,tol=1e-1,maxIter=10,out=-1) ->
                  blockBiCGSTB(A,B,M1=M,tol=tol,x=X,maxIter=maxIter,out=out)
sBiCG = getBlockIterativeSolver(IterMethod,out=1);
X, = solveLinearSystem(Ar,B,sBiCG);
@test vecnorm(Ar*X-B)/vecnorm(B) < sBiCG.tol

println("===  Test for nonsymmetric matrices ====");
n = 100
m = 4
sBiCG.PC = :jac
sBiCG.doClear=true
sBiCG.out=0
A = sprandn(n,n,5/n) + 10*speye(n)
B = randn(n,m)
X, = solveLinearSystem(A,B,sBiCG)
@test vecnorm(A*X - B)/vecnorm(B) < sBiCG.tol

sBiCG.isTranspose=true
Xt, = solveLinearSystem(A',B,sBiCG)
@test vecnorm(A*X - B)/vecnorm(B) < sBiCG.tol
@test vecnorm(Xt - X)/vecnorm(X) < 1e-13
