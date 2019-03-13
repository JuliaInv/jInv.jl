
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

pMisRefs    = Array{RemoteChannel}(undef, 6)
workerList  = workers()
nw          = nworkers()
pMisRefs[1] = initRemoteChannel(getMisfitParam,workerList[1%nw + 1],
                                pFor1,Wd1,bd1,SSDFun,fMod,gl1)
pMisRefs[2] = initRemoteChannel(getMisfitParam,workerList[2%nw + 1],
                                pFor2,Wd2,bd2,SSDFun,fMod,gl2)

pForRefs = Array{RemoteChannel}(undef, 4)
pForRefs[1] = initRemoteChannel(LSparam,workerList[3%nw+1],A[i3,:],[])
pForRefs[2] = initRemoteChannel(LSparam,workerList[4%nw+1],A[i4,:],[])
pForRefs[3] = initRemoteChannel(LSparam,workerList[5%nw+1],A[i3,:],[])
pForRefs[4] = initRemoteChannel(LSparam,workerList[6%nw+1],A[i4,:],[])

speye_t = (n)->SparseMatrixCSC(I, n, n);

Mesh2Mesh = Array{Future}(undef, 2)
for i=1:length(Mesh2Mesh)
	k = pForRefs[i].where
	if i==1
		Mesh2Mesh[i] = remotecall(speye_t,k,Minv.nc);
		wait(Mesh2Mesh[i]);
	elseif i==2
		Mesh2Mesh[i] = remotecall(identity,k,1.0);
		wait(Mesh2Mesh[i]);
	end
end

Wd = Array{Array{Float64}}(undef, 2);
dobs = Array{Array{Float64}}(undef, 2);
dobs[1] = bd3;
dobs[2] = bd4;
Wd[1] = Wd3;
Wd[2] = Wd4;
pMisRefs[3:4] = getMisfitParam(pForRefs[1:2],Wd,dobs,SSDFun,Iact,sigmaBack,Mesh2Mesh);
pMisRefs[5:6] = getMisfitParam(pForRefs[3:4],Wd,dobs,SSDFun,Iact,sigmaBack); # single mesh

M2M3		 = sparse(1.0I,Minv.nc,Minv.nc);
M2M4		 = sparse(1.0I,Minv.nc,Minv.nc);
pMis    = Array{MisfitParam}(undef, 6)
pMis[1] = getMisfitParam(pFor1,Wd1,bd1,SSDFun,fMod,gl1);
pMis[2] = getMisfitParam(pFor2,Wd2,bd2,SSDFun,fMod,gl2);
pMis[3] = fetch(pMisRefs[3]);
pMis[4] = fetch(pMisRefs[4]);
pMis[5] = fetch(pMisRefs[5]);
pMis[6] = fetch(pMisRefs[6]);

# setup pInv
alpha        = 1.
x0           = (sum(xtrue)./length(xtrue))*ones(Minv.nc)
boundsLow    = minimum(xtrue)*ones(Minv.nc)
boundsHigh   = maximum(xtrue)*ones(Minv.nc)
sigmaBack    = zeros(Minv.nc)

#  solve with automatic distribution
pInv         = getInverseParam(Minv,fMod,diffusionReg,alpha,x0,boundsLow,boundsHigh)
pInv.maxIter = 5
x1, = projGN(x0,pInv,pMis)

pInv.maxIter = 5
x2, = projGN(x0,pInv,pMisRefs)
@test norm(x1-x2)/norm(x1) .< 1e-12

clear!(pMisRefs);
for k=1:length(pMis)
	clear!(pMis[k],clearPFor=true, clearData=true,clearMesh2Mesh=true);
end

#Finalize pFor references
for pFo in pForRefs
	finalize(pFo)
end
