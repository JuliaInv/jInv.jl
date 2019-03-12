export getData


"""
function jInv.ForwardShare.getData(sigma, pFor, ...)

solves forward model and computes data. getData is application specific and some guidelines
how to create getData for a new problem can be found in examples/tutorialBuildYourOwn.ipynb

Inputs:

   sigma  -  current parameters
   pFor   -  description of forward problems (ForwardProbType, Array{ForwardProbType}, Array{Future}

   Some methods of getData require further arguments.

Output:

    Dobs   - simulated data
    pFor   - modified forward problem type

"""
function getData(sigma::Union{RemoteChannel, Vector},
                 pFor::ForwardProbType,
                 Mesh2Mesh::Mesh2MeshTypes,
                 doClear::Bool=false)
#=
    computations on one worker
=#
    sig = interpGlobalToLocal(fetch(sigma),fetch(Mesh2Mesh))
    dobs,pFor   = getData(sig,pFor,doClear)
    Dobs  = remotecall(identity,myid(),dobs)
    return Dobs,pFor
end

function getData(sigma::Union{RemoteChannel, Vector},
                 pFor::RemoteChannel,
                 Mesh2Mesh::Mesh2MeshTypes,
                 doClear::Bool=false)
#=
    modify pFor on current worker (e.g., to keep factorizations)
=#
    pF = take!(pFor)
    Dobs, pF = getData(fetch(sigma), pF, fetch(Mesh2Mesh), doClear)
    put!(pFor,pF)
    return Dobs,pFor
end

function getData(
                 sigma::Vector,
                 pFor::Array{FPT},
                 Mesh2Mesh::Array{T}=ones(length(pFor)),
                 doClear::Bool=false,
                 workerList::Vector=workers()) where {FPT<:ForwardProbType,T<:Mesh2MeshTypes}
#=
    parallel forward simulation with dynamic scheduling (i.e., elements in pFor get sent to remote workers on the fly)
=#
    i=1; nextidx() = (idx = i; i+=1; idx)

    Dobs = Array{Any}(undef,length(pFor))
    workerList = intersect(workers(),workerList)
    if isempty(workerList)
        error("getData: workers do not exist!")
    end
    @sync begin
        for p=workerList
            @async begin
                while true
                    idx = nextidx()
                    if idx > length(pFor); break; end
                    Dobs[idx],pFor[idx] = remotecall_fetch(getData,p,sigma,pFor[idx],Mesh2Mesh[idx],doClear)
                end
            end
        end
    end
    return Dobs,pFor
end

function getData(
                 sigma::Vector,
                 pFor::Array{RemoteChannel},
                 Mesh2Mesh::Array{T}=ones(length(pFor)),
                 doClear::Bool=false) where {T<:Mesh2MeshTypes}
#=
    parallel forward simulation with static scheduling (i.e., pFors are distributed a-priorily)
=#

    Dobs = Array{Future}(undef, length(pFor))
    workerList = getWorkerIds(pFor)
    sigmaRef = Array{RemoteChannel}(undef, maximum(workers()))
    @sync begin
        for p=workerList
            @async begin
                sigmaRef[p] = initRemoteChannel(identity,p,sigma)  # send model to worker
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
