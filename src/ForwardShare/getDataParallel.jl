export getData



function getData(sigma::Future,pFor::RemoteChannel,
                 Mesh2Mesh::Union{RemoteChannel, Future,SparseMatrixCSC,AbstractFloat},
                 doClear::Bool=false)
	#=
		load a forward problem from RemoteRef
	=#
	# get mesh 2 mesh interpolation matrix
	sig = interpGlobalToLocal(fetch(sigma),fetch(Mesh2Mesh))
	pF  = take!(pFor)
	dobs,pF   = getData(sig,pF,doClear)
	put!(pFor,pF)
	Dobs  = remotecall(identity,myid(),dobs)
	return Dobs,pFor
end

function getData(sigma::Vector,pFor::ForwardProbType,Mesh2Mesh::Union{SparseMatrixCSC,AbstractFloat},
                 doClear::Bool=false)

	# get mesh 2 mesh interpolation matrix
	sig = interpGlobalToLocal(sigma,Mesh2Mesh)
	dobs,pFor   = getData(sig,pFor,doClear)
	Dobs  = remotecall(identity,myid(),dobs)
	return Dobs,pFor
end

function getData{FPT<:ForwardProbType,T<:Union{Future,RemoteChannel,SparseMatrixCSC,AbstractFloat}}(
                 sigma::Vector,
                 pFor::Array{FPT},
                 Mesh2Mesh::Array{T}=ones(length(pFor)),
                 doClear::Bool=false,
                 workerList::Vector=workers())
	#=
		load and solve forward problems in parallel
	=#
	i=1; nextidx() = (idx = i; i+=1; idx)
	## Compute Dobs
	Dobs = Array(Any,length(pFor))
	workerList = intersect(workers(),workerList)
	if isempty(workerList)
		error("getData: workers do not exist!")
	end
	@sync begin
		for p=workerList
			@async begin
				while true
					idx = nextidx()
					if idx > length(pFor)
						break
					end
#					P = Mesh2Mesh[idx]
# 					if isa(P,SparseMatrixCSC) && eltype(P.nzval)==Int16
# 						P = SparseMatrixCSC(P.m,P.n,P.colptr,P.rowval,2.^float(3*P.nzval))
# 					end
					Dobs[idx],pFor[idx]    = remotecall_fetch(getData,p,sigma,pFor[idx],Mesh2Mesh[idx],doClear)
				end
			end
		end
	end
	return Dobs,pFor
end

function getData{T<:Union{RemoteChannel,Future,SparseMatrixCSC,AbstractFloat}}(
                 sigma::Vector,
                 pFor::Array{RemoteChannel},
                 Mesh2Mesh::Array{T}=ones(length(pFor)),
                 doClear::Bool=false)
#=
    parallel forward simulation with static scheduling (i.e., pFors are distributed a-priorily)
=#

    Dobs = Array(Future,length(pFor))
    workerList = getWorkerIds(pFor)
    sigmaRef = Array(Future,maximum(workers()))
    @sync begin
        for p=workerList
            @async begin
                sigmaRef[p] = remotecall(identity,p,sigma)  # send model to worker
                for idx=1:length(pFor)                      
                    if p==pFor[idx].where                  
                        Dobs[idx],pFor[idx] = remotecall_fetch(getData,p,sigmaRef[p],pFor[idx],Mesh2Mesh[idx],doClear)
                    end
                end
            end
        end
    end
    return Dobs,pFor
end
