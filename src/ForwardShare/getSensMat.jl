export getSensMat

"""
S = function getSensMat(...)

constructs sensitivity matrix.

WARNING: For large-scale problems this will be prohibively
         expensive. Use with caution

Inputs:

	m    - model
	pFor - forward problems

Examples:

	S = getSensMat(m,pFor)            # single pFor
	S = getSensMat(m,[pFor1;pFor2])   # multiple pFor's
	S = getSensMat(m,pForRef)         # pFor as remote reference
"""
function getSensMat(w::Vector,pFor::ForwardProbType)
	(n,m) = getSensMatSize(pFor)
	J    = zeros(n,m)

	if min(m,n) > 1e4; error("sensitivity matrix is too big to build column- or rowwise."); end

	if  n < m # decide which way is less work
		I = eye(n)
		for k=1:n
			tv = vec(getSensTMatVec(vec(I[:,k]),w,pFor))
			J[k,:] = tv
		end
	else
		I = eye(m)
		for k=1:m
			J[:,k] = vec(getSensMatVec(vec(I[:,k]),w,pFor))
		end
	end
	return J
end

function getSensMat(w,pFor::RemoteChannel)
	if pFor.where != myid()
		return remotecall_fetch(getSensMat,pFor.where,w,pFor)
	end
	return getSensMat(fetch(w),fetch(pFor))
end

getSensMat(w::Future,pFor::ForwardProbType) = getSensMat(fetch(w),pFor)

function getSensMat(m,pFor::Array,workerList=workers())

	S = Array{Any}(length(pFor))
	i=1; nextidx() = (idx = i; i+=1; idx)

	workerList = intersect(workers(),workerList)
	if isempty(workerList)
		error("getSensMat: specified workers do not exist!")
	end

	mRef = Array{Future}(maximum(workers()))

	@sync begin
		for p=workerList
			@async begin
				mRef[p] = remotecall(identity,p,m)
				while true
					idx = nextidx()
					if idx > length(pFor)
						break
					end
					S[idx]    = remotecall_fetch(getSensMat,p,mRef[p],pFor[idx])
				end
			end
		end
	end
	return vcat(tuple(S...)...)
end
