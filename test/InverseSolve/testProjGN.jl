include("../setupTests.jl")
# build domain and true image
domain = [0.0 1.0 0.0 1.0]
n      = [16,16]
Minv   = getRegularMesh(domain,n)
xc     = getCellCenteredGrid(Minv)
xtrue = sin(2*pi*xc[:,1]).*sin(pi*xc[:,2])

# get noisy data 
A     = sparse(1.0I,Minv.nc,Minv.nc)
ids   = sort(rand(1:Minv.nc,round(Int64,Minv.nc*.8)))
A     = A[ids,:]
btrue = A*xtrue
bdata = btrue + .1*randn(size(btrue))/norm(btrue,Inf)
Wd    = ones(length(btrue))

# generate misfit param
i1     = (1:round(Int64,size(A,1)/2))
i2     = (round(Int64,size(A,1)/2)+1:size(A,1))
pFor1  = LSparam(A[i1,:],[])
pFor2  = LSparam(A[i2,:],[])
sigmaBack = zeros(length(xtrue))
Iact         = sparse(1.0I,Minv.nc,Minv.nc)
gl1          = getGlobalToLocal(1.0,sigmaBack)
gl2          = getGlobalToLocal(Iact,sigmaBack)
bd1          = bdata[i1]
bd2          = bdata[i2]
Wd1          = Wd[i1]
Wd2          = Wd[i2]

pMisRefs    = Array(RemoteChannel,2)
pMisRefs[1] = initRemoteChannel(getMisfitParam,workers()[1],pFor1,Wd1,bd1,SSDFun,fMod,gl1) 
pMisRefs[2] = initRemoteChannel(getMisfitParam,workers()[min(2,nworkers())],pFor2,Wd2,bd2,SSDFun,fMod,gl2) 

pMis    = Array{MisfitParam}(2)
pMis[1] = getMisfitParam(pFor1,Wd1,bd1,SSDFun,fMod,gl1) 
pMis[2] = getMisfitParam(pFor2,Wd2,bd2,SSDFun,fMod,gl2) 

# setup pInv
alpha        = 1.
x0           = mean(xtrue)*ones(Minv.nc)
boundsLow    = -.5*ones(Minv.nc)
boundsHigh   = 0.5*ones(Minv.nc)
sigmaBack    = zeros(Minv.nc)

print("\t== compare with projGN...")
pInv         = getInverseParam(Minv,fMod,diffusionReg,alpha,x0,boundsLow,boundsHigh)
pInv.maxIter = 5
x1t, = projGN(copy(x0),pInv,pMis)
pInv.maxIter = 5
x2t, = projGN(copy(x0),pInv,pMisRefs,out=1)
@test norm(x1t-x2t)/norm(x1t) < 1e-12
@test all(x1t.>=boundsLow)
@test all(x2t.>=boundsLow)
@test all(x1t.<=boundsHigh)
@test all(x2t.<=boundsHigh)
print("passed! ===")

print("\t== test projGN with projCG solver...")
pInv.maxIter = 5
x1, = projGN(copy(x0),pInv,pMis)
pInv.maxIter = 5
x2, = projGN(copy(x0),pInv,pMisRefs,out=1)
@test norm(x1-x2)/norm(x1) < 1e-12
@test norm(x1t-x1)/norm(x1t) < 1e-12
@test norm(x2t-x2)/norm(x2t) < 1e-12
@test all(x1.>=boundsLow)
@test all(x2.>=boundsLow)
@test all(x1.<=boundsHigh)
@test all(x2.<=boundsHigh)
print("passed! ===")

print("\t== test projGN with direct solver...")
pInv.maxIter = 5
x1, = projGN(copy(x0),pInv,pMis,solveGN=projGNexplicit)
pInv.maxIter = 5
x2, = projGN(copy(x0),pInv,pMisRefs,solveGN=projGNexplicit,out=1)
@test norm(x1-x2)/norm(x1) < 1e-12
@test norm(x1t-x1)/norm(x1t) < 1e-2
@test norm(x2t-x2)/norm(x2t) < 1e-2
@test all(x1.>=boundsLow)
@test all(x2.>=boundsLow)
@test all(x1.<=boundsHigh)
@test all(x2.<=boundsHigh)

print("\t== test projGN with singular system...")
pInv.alpha = 0.
x2, = projGN(copy(x0),pInv,pMis[1],solveGN=projGNexplicit,out=1)
@test all(x2.>=boundsLow)
@test all(x2.<=boundsHigh)


print("passed! ===")

