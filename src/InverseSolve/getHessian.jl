using jInv.InverseSolve
export getHessian

"""
function H = getHessian(...)

builds and returns Hessian matrix of misfit functions.

Input:

	sig   - model
	pMis  - misfit param (supported: single, array of misfit params, remote references)
	d2F   - Hessian of distance (supported: single, array, remote references)

Output:

	H     - Hessian

Examples:

	H  = getHessian(sigma, pMis, d2F)
	Ha = getHessian(sigma, [pMis1;pMis2], [d2F1; d2F2])

"""
function getHessian(sig::Vector,  # conductivity on inv mesh (active cells only)
                    pMis::MisfitParam,
                    d2F)
	sigma,dsigma = pMis.modelfun(sig)

	P = pMis.gloc.PForInv
	J = getSensMat(sigma,pMis.pFor)

	dr = J*P*dsigma

	return dr'*sdiag(d2F)*dr
end

function getHessian(sig,pMis::Array,d2F::Array,workerList=workers())

	Hs = Array{Any}(undef, length(pMis))
	i=1; nextidx() = (idx = i; i+=1; idx)

	workerList = intersect(workers(),workerList)
	if isempty(workerList)
		error("getHessian: specified workers do not exist!")
	end

	sigRef = Array{Future}(undef, maximum(workerList))

	@sync begin
		for p=workerList
			@async begin
				sigRef[p] = remotecall(identity,p,sig)
				while true
					idx = nextidx()
					if idx > length(pMis)
						break
					end
					Hs[idx]    = remotecall_fetch(getHessian,p,sigRef[p],pMis[idx],d2F[idx])
				end
			end
		end
	end
	H = Hs[1]
	for k=2:length(Hs)
		H += Hs[k]
	end

	return H
end

getHessian(sig::Union{RemoteChannel,Future}, pMis::Union{RemoteChannel,Future}, d2F::Union{RemoteChannel,Future}) = getHessian(fetch(fetch(sig)),fetch(pMis),fetch(d2F))
getHessian(sig::Array, pMis::Union{RemoteChannel,Future}, d2F::Union{RemoteChannel,Future}) = getHessian(sig,fetch(pMis),fetch(d2F))
getHessian(sig::Union{RemoteChannel,Future}, pMis::MisfitParam, d2F::Array) = getHessian(fetch(sig), pMis, d2F)
