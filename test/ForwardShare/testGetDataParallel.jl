
# build domain and true image
domain = [0.0 1.0 0.0 1.0]
n      = [16,16]
Minv   = getRegularMesh(domain,n)
xc     = getCellCenteredGrid(Minv)
xtrue  = sin.(2*pi*xc[:,1]).*sin.(pi*xc[:,2])
A      = sparse(1.0I,Minv.nc,Minv.nc)
ids    = sort(rand(1:Minv.nc,round(Int64,Minv.nc*.8)))
A      = A[ids,:]
btrue  = A*xtrue

#Generate two forward problems
i1     = (1:round(Int64,size(A,1)/2))
i2     = (round(Int64,size(A,1)/2)+1:size(A,1))
pFor1  = LSparam(A[i1,:],[])
pFor2  = LSparam(A[i2,:],[])
sigmaBack    = zeros(length(xtrue))
Iact         = sparse(1.0I,Minv.nc,Minv.nc)
gl1          = getGlobalToLocal(1.0,sigmaBack)
gl2          = getGlobalToLocal(Iact,sigmaBack)

#Test with locally stored forward problems
localPfors = [pFor1;pFor2]
D1, = getData(xtrue,localPfors)
@test norm([fetch(D1[1]);fetch(D1[2])]-btrue)/norm(btrue) < 1e-10

#Test with forward problems on remote workers
pForRefs    = Array{RemoteChannel}(undef,2)
workerList  = workers()
nw          = nworkers()
pForRefs[1] = initRemoteChannel(LSparam,workerList[1%nw+1],A[i1,:],[])
pForRefs[2] = initRemoteChannel(LSparam,workerList[2%nw+1],A[i2,:],[])
D2,pForRefs = getData(xtrue,pForRefs,ones(length(pForRefs)),true)
@test norm([fetch(D2[1]);fetch(D2[2])]-btrue)/norm(btrue) < 1e-10

#Test doClear flag.
@test fetch(pForRefs[1]).A    == sparse(I,0,0)
@test fetch(pForRefs[1]).Ainv == []
@test fetch(pForRefs[2]).A    == sparse(I,0,0)
@test fetch(pForRefs[2]).Ainv == []
@test_throws Exception D3, = getData(xtrue,pForRefs)
