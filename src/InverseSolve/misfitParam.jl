export MisfitParam,getMisfitParam

const ArrayUnionLocRemote = Union{Array{RemoteChannel},Array{Future},Array{Float64}}
const UnionLocRemote      = Union{RemoteChannel,Future,AbstractFloat}

"""
mutable struct jInv.InverseSolve.MisfitParam

Type storing information about one term in the misfit

F(m) = sum_i^n phi_i(pFor(model(m)),dobs,Wd)

Fields:

    pFor::ForwardProbType  - forward problem
    Wd                     - inverse standard deviation
    dobs                   - observed data
    misfit::Function       - misfit function
    modelfun::Function     - model function (evaluated locally)
    gloc                   - mapping from inverse to forward mesh (including active cells projection).

Constructors:

getMisfitParam(pFor,Wd,dobs,misfit,model,gloc=getGlobalToLocal(1.0))

getMisfitParam(pForRFs::Array{RemoteChannel}, Wd::Array, dobs::Array, misfit::Function,
                            Iact,sigmaBack::Vector,
                            Mesh2MeshRFs::Union{Array{RemoteChannel},Array{Float64}}=ones(length(pForRFs)),
                            modelfun::Function=identityMod,fname="")

"""
mutable struct  MisfitParam
    pFor::ForwardProbType
    Wd
    dobs
    misfit::Function
    modelfun::Function
    gloc
end

"""
function jInv.InverseSolve.getMisfitParam

Required Input:

Input for a call:
    getMisfitParam(pFor,Wd,dobs,misfit,model,gloc=getGlobalToLocal(1.0))

    pFor::ForwardProbType  - forward problem
    Wd                     - inverse standard deviation of data
    dobs                   - observed data
    misfit::Function       - Misfit function
    modelfun::Function     - model function. If all misfits have same model function,
                             put modelfun in inverseParam and identityMod here for efficiency.
    gloc                   - mapping from inverse to forward mesh (default: identity). See globalToLocal.jl

Input for a call:
    getMisfitParam(pForRFs::Array{RemoteChannel}, Wd::Array, dobs::Array, misfit::Function,
                            Iact,sigmaBack::Vector, Mesh2MeshRFs::Union{Array{RemoteChannel},Array{Float64}}=ones(length(pForRFs)),
                            modelfun::Function=identityMod,fname="")

    pForRFs::Array{RemoteChannel}                      - Array of remote references for forward problem parameters on each worker.
    Wd::Array                                - Array of inverse standard deviation for the data on each worker.
    dobs::Array                                - Array of observed data for each worker.
    misfit::Function                            - Misfit function
    Iact                                        - Projector to active cells.
    sigmaBack::Vector                            - Background model ("frozen" cells).
    Mesh2MeshRFs::Union{Array{RemoteChannel},Array{Float64}}    - mapping from inverse to forward mesh (default: identity)
    modelfun::Function=identityMod                         - model function. If all misfits have same model function,
                                          put modelfun in inverseParam and identityMod here for efficiency.
    fname=""                                 - (optional)

"""
function getMisfitParam(pFor::ForwardProbType, Wd, dobs, misfit::Function, modelfun::Function=identityMod, gloc=getGlobalToLocal(1.0))
    return MisfitParam(pFor,Wd,dobs,misfit,modelfun,gloc);
end

function getMisfitParam(pForRFs::Array{RemoteChannel}, Wd::Array, dobs::Array, misfit::Function,
                            Iact,sigmaBack::Union{Vector,AbstractFloat,AbstractModel},
                            Mesh2MeshRFs::ArrayUnionLocRemote=ones(length(pForRFs)),
                            modelfun::Function=identityMod,fname="")
    pMis                = Array{RemoteChannel}(undef, length(pForRFs));
    # the next references are for not to broadcast Iact and sigmaBackground
    IactRFs            = Array{Future}(undef,maximum(workers()));
    sigmaBackgroundRFs = Array{Future}(undef,maximum(workers()));

    @sync begin
        for p=workers()
            @async begin
                # communicate Iact and sigmaBack
                IactRFs[p]            = remotecall_wait(identity,p,Iact)            # send active cells projection matrix
                sigmaBackgroundRFs[p] = remotecall_wait(identity,p,sigmaBack) # send background conductivity
                for idx=1:length(pForRFs)
                    if p==pForRFs[idx].where
                        pMis[idx] = initRemoteChannel(getMisfitParam,p,pForRFs[idx],Wd[idx],dobs[idx],misfit,IactRFs[p],sigmaBackgroundRFs[p],Mesh2MeshRFs[idx],modelfun,fname);
                    end
                end
            end
        end
    end
    return pMis;
end

function getMisfitParam(pForRF::RemoteChannel, Wd, dobs, misfit::Function,
                            Iact::Future,sigmaBack::Future,
                            Mesh2MeshRF::UnionLocRemote = 1.0,
                            modelfun::Function=identityMod,fname="")
    # Single worker version of getMisfitParam.
    worker = pForRF.where;
    if !isa(Mesh2MeshRF,AbstractFloat)
        Mesh2Mesh     = fetch(Mesh2MeshRF);
        finalize(Mesh2MeshRF)  # to prevent memory leak
        if worker!=Mesh2MeshRF.where
            error("getMisfitParam::Mesh2Mesh and pFor not on the same worker");
        end
    else
        Mesh2Mesh = Mesh2MeshRF;
    end
    pFor         = fetch(pForRF);
    Iact         = fetch(Iact);
    sigmaBack     = fetch(sigmaBack);
    Wd              = fetch(Wd);
    dobs            = fetch(dobs);
    pMis         = getMisfitParam(pFor,Wd,dobs,misfit,modelfun, prepareGlobalToLocal(Mesh2Mesh,Iact,sigmaBack,fname));
    return pMis;
end

import jInv.Utils.clear!
function clear!(P::MisfitParam;clearPFor::Bool=true, clearData::Bool=false,clearMesh2Mesh::Bool=false)
    if clearPFor
        clear!(P.pFor);
    end
    if clearData
        P.Wd   = 0;
        P.dobs = 0;
    end
    if clearMesh2Mesh
        P.gloc = getGlobalToLocal(1.0);
    end
    return P;
end
