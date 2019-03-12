export prepareMesh2Mesh

"""
function jInv.ForwardShare.prepareMesh2Mesh(Mloc,Minv,compact)

builds interpolation matrix from inverse mesh, Minv, to PDE mesh, Mloc.
Used for efficiency in mesh-decoupling.

Examples:

    P = prepareMesh2Mesh(Msmall, Mbig, true)
    modelLoc =  interpGlobalToLocal(modelInv,P) # uncompresses the indices and does P'*modelInv
    kwargs = Additional arguments used by getInterpolationMatrix

Inputs:

    Mloc    - pde mesh (can also be remote ref or part of pFor)
    Minv    - inverse mesh
    compact - flag to compress indices/values of interpolation matrix (default: true)

Outputs:

    P'       - interpolation matrix (sparse). To reduce storage, transpose is returned.

"""
function prepareMesh2Mesh(Mfwd::AbstractMesh, Minv::AbstractMesh, compact::Bool=true; kwargs...)

    P = sparse(getInterpolationMatrix(Minv,Mfwd; kwargs...)');
    if compact
        Ps = SparseMatrixCSC(P.m,P.n,round.(UInt32,P.colptr),round.(UInt32,P.rowval),round.(Int8,log2.(P.nzval)/3))
    else
        Ps = SparseMatrixCSC(P.m,P.n,round.(UInt32,P.colptr),round.(UInt32,P.rowval),P.nzval)
    end
    return Ps
end

function prepareMesh2Mesh(pF::ForwardProbType, Minv::AbstractMesh, compact::Bool=true; kwargs...)
    return prepareMesh2Mesh(pF.Mesh, Minv, compact; kwargs...)
end

function prepareMesh2Mesh(pFor::RemoteChannel, MinvRef::Future, compact::Bool=true; kwargs...)
    return prepareMesh2Mesh(fetch(pFor), fetch(MinvRef), compact; kwargs...)
end

function prepareMesh2Mesh(pFor::Array{RemoteChannel},Minv::AbstractMesh,compact::Bool=true; kwargs...)

    Mesh2Mesh = Array{RemoteChannel}(undef, length(pFor))
    # find out which workers are involved
    workerList = []
    for k=1:length(pFor)
        push!(workerList,pFor[k].where)
    end
    workerList = unique(workerList)
    # send sigma to all workers
    MinvRef = Array{Future}(undef, maximum(workers()))

    @sync begin
        for p=workerList
            @async begin
                MinvRef[p] = remotecall_wait(identity,p,Minv)   # send model to workers
                for idx=1:length(pFor)
                    if p==pFor[idx].where
                        Mesh2Mesh[idx] = initRemoteChannel(prepareMesh2Mesh,p,pFor[idx],MinvRef[p],compact; kwargs...)
                    end
                end
            end
        end
    end
    return Mesh2Mesh
end
