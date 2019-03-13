
# get random matrix for least squares problem
A     = ( sprandn(100,13,.1), sprandn(100,200,.2), randn(100,24))
i1     = (1:round(Int64,size(A[1],1)/2))
i2     = (round(Int64,size(A[1],1)/2)+1:size(A[1],1))
pFor1  = LSparam(A[1][i1,:],[])
pFor2  = LSparam(A[1][i2,:],[])

# single LS problem
pFor   = LSparam(A[2],[])
# parallel with dynamic scheduling
pForp  = [pFor1; pFor2];
# remote references
pForRef    = Array{RemoteChannel}(undef,2)
i1     = (1:round(Int64,size(A[3],1)/2))
i2     = (round(Int64,size(A[3],1)/2)+1:size(A[3],1))
A1     = A[3][i1,:]
A2     = A[3][i2,:]
workerList = workers()
nw         = length(workerList)
pForRef[1] = initRemoteChannel(LSparam,workerList[1%nw+1],A1,[])
pForRef[2] = initRemoteChannel(LSparam,workerList[2%nw+1],A2,[])

pFors = (pForp,pFor,pForRef)

for k=1:length(pFors)
	# test getSensMatVec for pFor as $(typeof(pFors[k]))...
	(mt,nt) = getSensMatSize(pFors[k])
	nd      = getNumberOfData(pFors[k])
	(mk,nk) = size(A[k])

	@test mt==nd
	@test mk==mt
	@test nt==nk

	At = getSensMat(randn(nk),pFors[k])
	if At[1] isa Future
	    At = vcat([fetch(Ai) for Ai in At]...)
	end
	@test norm(At-A[k],Inf)/norm(A[k],Inf) < 1e-14
end
