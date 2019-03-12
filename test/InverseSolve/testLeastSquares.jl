
# build domain and true image
domain = [0.0 1.0 0.0 1.0]
n      = [16,16]
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

pMisRefs    = Array{RemoteChannel}(undef, 2)
workerList  = workers()
nw          = nworkers()
pMisRefs[1] = initRemoteChannel(getMisfitParam,workerList[1%nw+1],pFor1,Wd1,bd1,SSDFun,fMod,gl1)
pMisRefs[2] = initRemoteChannel(getMisfitParam,workerList[2%nw+1],pFor2,Wd2,bd2,SSDFun,fMod,gl2)

pMis    = Array{MisfitParam}(undef, 2)
pMis[1] = getMisfitParam(pFor1,Wd1,bd1,SSDFun,fMod,gl1)
pMis[2] = getMisfitParam(pFor2,Wd2,bd2,SSDFun,fMod,gl2)

# setup pInv
alpha        = 1.
x0           = (sum(xtrue)/length(xtrue))*ones(Minv.nc)
boundsLow    = minimum(xtrue)*ones(Minv.nc)
boundsHigh   = maximum(xtrue)*ones(Minv.nc)
sigmaBack    = zeros(Minv.nc)

# test projGN
pInv         = getInverseParam(Minv,fMod,diffusionReg,alpha,x0,boundsLow,boundsHigh)
pInv.maxIter = 5
x1, = projGN(x0,pInv,pMis)
pInv.maxIter = 5
x2, = projGN(x0,pInv,pMisRefs,out=1)
@test norm(x1-x2)/norm(x1) < 1e-12
@test all(x1.>=boundsLow)
@test all(x2.>=boundsLow)
@test all(x1.<=boundsHigh)
@test all(x2.<=boundsHigh)
# test barrierGNCG...")

x1, = barrierGNCG(x0,pInv,pMis);
x2, = barrierGNCG(x0,pInv,pMisRefs,out=1);
@test norm(x1-x2)/norm(x1) < 1e-12
@test all(x1.>=boundsLow)
@test all(x2.>=boundsLow)
@test all(x1.<=boundsHigh)
@test all(x2.<=boundsHigh)


# test iteratedTikhonov...")
pInv.maxIter = 2
pInv.alpha   = 100.
nAlpha = 3
alphaFac = 10.
alphaMin = alpha/(alphaFac^nAlpha)
alphaParam = [alpha;alphaMin;alphaFac;nAlpha]
targetMisfit = 20.0
x1,Dc,flag1,     = iteratedTikhonov(x0,pInv,pMis,alphaParam,targetMisfit)

pInv.alpha = 100.
pInv.mref  = x0
x2,Dc,flag2,hist = iteratedTikhonov(x0,pInv,pMisRefs,alphaParam,targetMisfit,solveGN=projGNexplicit)
@test typeof(hist) <: Array{InverseSolve.GNhis}
@test norm(x1-x2)/norm(x1) < 1e-12
@test all(x1.>=boundsLow)
@test all(x2.>=boundsLow)
@test all(x1.<=boundsHigh)
@test all(x2.<=boundsHigh)

# test projSD...")
pInv.maxIter = 20
x1, = projSD(x0,pInv,pMis)
pInv.maxIter = 20
x2, = projSD(x0,pInv,pMisRefs,out=1)
@test norm(x1-x2)/norm(x1) < 1e-12
@test all(x1.>=boundsLow)
@test all(x2.>=boundsLow)
@test all(x1.<=boundsHigh)
@test all(x2.<=boundsHigh)
