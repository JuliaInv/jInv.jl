export getSensMat

"""
S = function getSensMat(...)

constructs sensitivity matrix.

WARNING: For large-scale problems this will be prohibively
         expensive. Use with caution

Inputs:

	sigma - model
	pFor  - forward problems

Examples:

	S = getSensMat(sigma, pFor)            # single pFor

    Some methods of getData require further arguments.

"""
function getSensMat(sigma::Vector, 
					pFor::ForwardProbType)

	(n, m) = getSensMatSize(pFor)
	sensMat = zeros(n, m)

	if min(m, n) > 1e4
		error("sensitivity matrix is too big to build column- or rowwise.")
	end

	if  n < m # decide which way is less work
		I = eye(n)
		for k = 1:n
			tv = vec(getSensTMatVec(vec(I[:, k]), sigma, pFor))
			sensMat[k, :] = tv
		end
	else
		I = eye(m)
		for k = 1:m
			sensMat[:, k] = vec(getSensMatVec(vec(I[:, k]), sigma, pFor))
		end
	end
	return sensMat
end

function getSensMat(sigma::Union{RemoteChannel, Vector},
                    pFor::ForwardProbType,
                    Mesh2Mesh::Mesh2MeshTypes)

    sig = interpGlobalToLocal(fetch(sigma), fetch(Mesh2Mesh))
    sensMat = getSensMat(sig, pFor)
    sensMat = remotecall(identity, myid(), sensMat)
    return sensMat
end

function getSensMat(sigma::Union{RemoteChannel, Vector},
                    pFor::RemoteChannel,
                    Mesh2Mesh::Mesh2MeshTypes=1.)

    pF = take!(pFor)
    sensMat = getSensMat(fetch(sigma), pF, fetch(Mesh2Mesh))
    put!(pFor, pF)
    return sensMat
end

function getSensMat{T<:Mesh2MeshTypes}(sigma::Vector,
                 					   param::Array{RemoteChannel},
                 					   Mesh2Mesh::Array{T}=ones(length(param)))

    sensMat = Array{Future}(length(param))
    workerList = getWorkerIds(param)
    sigmaRef = Array{RemoteChannel}(maximum(workers()))
    @sync begin
        for p=workerList
            @async begin
                sigmaRef[p] = initRemoteChannel(identity, p, sigma)  # send model to worker
                for idx = 1:length(param)
                    if p == param[idx].where
                        sensMat[idx] = remotecall_fetch(getSensMat, p, sigmaRef[p], param[idx], Mesh2Mesh[idx])
                    end
                end
            end
        end
    end
    return sensMat
end
