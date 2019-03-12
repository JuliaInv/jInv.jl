using Test
using LinearAlgebra
using SparseArrays
using Distributed

using jInv.Utils
using jInv.Mesh
using jInv.InverseSolve

# build domain and true image
domain = [0.0 1.0 0.0 1.0]
n      = [32,32]
Minv   = getRegularMesh(domain,n)
xc     = getCellCenteredGrid(Minv)
xtrue = sin.(2*pi*xc[:,1]).*sin.(pi*xc[:,2])

# get noisy data
A     = sparse(1.0I,Minv.nc,Minv.nc)
ids   = sort(rand(1:Minv.nc,round(Int64,Minv.nc*.8)))
A     = A[ids,:]
btrue = A*xtrue
bdata = btrue + .1*randn(size(btrue))/norm(btrue,Inf)
Wd    = ones(length(btrue))

# generate misfit param
i1     = (1:round(Int64,size(A,1)/4));
i2     = (round(Int64,size(A,1)/4)+1:round(Int64,2*size(A,1)/4));
i3 	   = (round(Int64,2*size(A,1)/4)+1 : round(Int64,3*size(A,1)/4));
i4 	   = (round(Int64,3*size(A,1)/4)+1 : size(A,1));
pFor1  = LSparam(A[i1,:],[])
pFor2  = LSparam(A[i2,:],[])



sigmaBack = zeros(length(xtrue))
Iact         = sparse(1.0I,Minv.nc,Minv.nc);
gl1          = getGlobalToLocal(1.0,sigmaBack)
gl2          = getGlobalToLocal(Iact,sigmaBack)

bd1          = bdata[i1]
bd2          = bdata[i2]
bd3          = bdata[i3]
bd4          = bdata[i4]
Wd1          = Wd[i1]
Wd2          = Wd[i2]
Wd3          = Wd[i3]
Wd4          = Wd[i4]

pMis =  getMisfitParam(pFor1,Wd1,bd1,SSDFun,fMod,gl1)

pMisRefs    = Array{RemoteChannel}(undef, 2)
pMisRefs[1] = initRemoteChannel(getMisfitParam,workers()[1],pFor1,Wd1,bd1,SSDFun,fMod,gl1)
pMisRefs[2] = initRemoteChannel(getMisfitParam,workers()[min(2,nworkers())],pFor2,Wd2,bd2,SSDFun,fMod,gl2)


pForRefs = Array{RemoteChannel}(undef, 4)
pForRefs[1] = initRemoteChannel(LSparam,workers()[1],A[i3,:],[])
pForRefs[2] = initRemoteChannel(LSparam,workers()[min(2,nworkers())],A[i4,:],[])
pForRefs[3] = initRemoteChannel(LSparam,workers()[min(3,nworkers())],A[i3,:],[])
pForRefs[4] = initRemoteChannel(LSparam,workers()[min(4,nworkers())],A[i4,:],[])


alpha        = 1.
x0           = (sum(xtrue)/length(xtrue))*ones(Minv.nc)

boundsLow    = minimum(xtrue)*ones(Minv.nc)
boundsHigh   = maximum(xtrue)*ones(Minv.nc)
pInv         = getInverseParam(Minv,fMod,diffusionReg,alpha,x0,boundsLow,boundsHigh)
display(pInv)

# Test with local misfit
# evaluate distance
sig,dsig = pInv.modelfun(x0)
Dc,F,dF,d2F,pMis,tMis = computeMisfit(sig,pMis,true)
# get Hessian
Hm = getHessian(sig,pMis,d2F)
# multiply with random direction and compare results
x = randn(size(Hm,2))
t1 = HessMatVec(x,pMis,sig,d2F)
t2 = Hm*x
@test norm(t1-t2,Inf)/norm(t1,Inf) < 1e-12

# Test with single remote ref misfit param
sig,dsig = pInv.modelfun(x0)
Dc,F,dF,d2F,pMisR,tMis = computeMisfit(sig,pMisRefs[1:1],true)
# get Hessian
Hm = getHessian(sig,pMisR[1],d2F[1])
# multiply with random direction and compare results
x = randn(size(Hm,2))
t1 = HessMatVec(x,pMisR,sig,d2F)
t2 = Hm*x
@test norm(t1-t2,Inf)/norm(t1,Inf) < 1e-12

# Test array of remote ref misfit param
# get Hessian
Hm = getHessian(sig,pMisR,d2F)
# multiply with random direction and compare results
x = randn(size(Hm,2))
t1 = HessMatVec(x,pMisR,sig,d2F)
t2 = Hm*x
@test norm(t1-t2,Inf)/norm(t1,Inf) < 1e-12

# Test array of misfit params
pMisA = [pMis;pMis;pMis]
sig,dsig = pInv.modelfun(x0)
Dc,F,dF,d2F,pMisA,tMis = computeMisfit(sig,pMisA,true)

# get Hessian
Hm = getHessian(sig,pMisA,d2F)
# multiply with random direction and compare results
x = randn(size(Hm,2))
t1 = HessMatVec(x,pMisA,sig,d2F)
t2 = Hm*x
@test norm(t1-t2,Inf)/norm(t1,Inf) < 1e-12
