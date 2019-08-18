
export HessMatVec

"""
function jInv.InverseSolve.HessMatVec

computes matrix-vector products with Hessian

    H(m)*x = J(m)*d2F*J(m)'*x

HessMatVec has different methods to make this efficient in the respective
cases. Type methods(HessMatVec) for a list.

"""
function HessMatVec(d2F::Union{Array{Float64,1},Array{Float32,1}}, x::Union{Array{Float64,1},Array{Float32,1}})
#=
   Hessian is real diagonal and represented as vector
=#
    return d2F .* x
end

function HessMatVec(d2F::Array{ComplexF64,1}, x::Array{ComplexF64,1})
#=
   Hessian is complex diagonal and represented as vector
=#
    return complex.(real.(d2F) .* real.(x), imag.(d2F) .* imag.(x))
end

function HessMatVec(d2F::SparseMatrixCSC{Float64}, x::Array{Float64,1})
#=
   Hessian is real and sparse
=#
    return d2F * x
end

function HessMatVec(d2F::SparseMatrixCSC{ComplexF64}, x::Array{ComplexF64,1})
#=
   Hessian is complex and sparse, vector is complex
=#
    return complex.(real.(d2F) * real.(x), imag.(d2F) * imag.(x))
end


function HessMatVec(d2F::SparseMatrixCSC{Float64}, x::Array{ComplexF64,1})
#=
   Hessian is real and sparse, vector is complex
=#
   return real2complex( d2F * complex2real(x))
end



function HessMatVec(x,
                    pMis::MisfitParam,
                    sig,  # conductivity on inv mesh (active cells only)
                    d2F)
#=
   Hessian includes forward problem all stored on current worker
=#
    # try
        sigma,dsigma = pMis.modelfun(sig)

        sigmaloc = interpGlobalToLocal(sigma,pMis.gloc.PForInv,pMis.gloc.sigmaBackground)
        xloc     = interpGlobalToLocal(dsigma*x,pMis.gloc.PForInv)
        Jx       = vec(getSensMatVec(xloc,sigmaloc,pMis.pFor))
        Jx       = HessMatVec(d2F,Jx)
        JTxloc   = getSensTMatVec(Jx,sigmaloc,pMis.pFor)
        JTx      = dsigma'*interpLocalToGlobal(JTxloc,pMis.gloc.PForInv) # =
        return JTx
    # catch err
        # if isa(err,InterruptException)
            # return -1
        # else
            # throw(err)
        # end
    # end
end


function HessMatVec(xRef::RemoteChannel,
                    pMisRef::RemoteChannel,
                    sigmaRef::RemoteChannel,
                    d2FRef::Future,
                    mvRef::RemoteChannel)
#=
   Hessian includes forward problem all stored as RemoteChannels or Futures
=#

    rrlocs = [xRef.where pMisRef.where sigmaRef.where d2FRef.where mvRef.where]
    if !all(rrlocs .== myid())
        println("WARNING: HessMatVec: Problem on worker ",myid()," not all remote refs are stored here, but rrlocs=",rrlocs);
    end

    # fetching and taking: should be no-ops as all RemoteRefs live on my worker
    t = time_ns();
    x     = fetch(xRef)
    sigma = fetch(sigmaRef)

    pMis  = take!(pMisRef)
    d2F   = fetch(d2FRef)
    commTime = (time_ns() - t)/1e+9;

    # compute HessMatVec

    compTime =@elapsed mvi  = HessMatVec(x,pMis,sigma,d2F)


    # putting: should be no ops.
    t = time_ns();
    mv   = take!(mvRef)
    put!(mvRef,mv+mvi)
    put!(pMisRef,pMis)
    commTime += (time_ns() - t)/1e+9;

    return true,commTime,compTime
end


function HessMatVec(x,
                    pMisRefs::Array{RemoteChannel,1},
                    sigma,
                    d2F,
                    indFors=1:length(pMisRefs))
#=
   Hessian includes multiple forward problem all stored as RemoteChannels or Futures
=#

    # find out which workers are involved
    workerList = getWorkerIds(pMisRefs)

    sigmaRef = Array{RemoteChannel}(undef,maximum(workers()))
    yRef = Array{RemoteChannel}(undef,maximum(workers()))
    zRef = Array{RemoteChannel}(undef,maximum(workers()))
    z = zeros(length(x))
    commTime = 0.0
    compTime = 0.0
    updateTimes(c1,c2) = (commTime+=c1; compTime+=c2)
    updateMV(x) = (z+=x)

    @sync begin
        for p = workerList
            @async begin
                # send model and vector to worker
                t = time_ns(); 
                sigmaRef[p] = initRemoteChannel(identity,p,sigma)
                yRef[p]     = initRemoteChannel(identity,p,x)
                zRef[p]     = initRemoteChannel(identity,p,zeros(length(x)))
                c1 = (time_ns() - t)/1e+9;
                updateTimes(c1,0.0)

                # do the actual computation
                for idx=indFors
                    if pMisRefs[idx].where==p
                        b,c1,c2 = remotecall_fetch(HessMatVec,p,yRef[p],pMisRefs[idx],sigmaRef[p],d2F[idx],zRef[p])
                        updateTimes(c1,c2)
                    end
                end

                # get back result (timing is a mixture of addition and computation)
                c1 = @elapsed updateMV(fetch(zRef[p]))
                updateTimes(c1,0.0)
            end
        end
    end
    matvc = z
    return matvc
end


function HessMatVec(x,pMis::Array{MisfitParam},sigma,d2F,indFors=1:length(pMis))
#=
   Hessian includes multiple forward problem all stored as RemoteChannels or Futures
=#
    numFor = length(pMis);

    # get process map
    i      = 1; nextidx = p -> i;
    procMap = zeros(Int64,numFor)
	
	###
	## Eran Treister:
	## The code below is a hack in case we're using a MUMPSsolver, who is no longer an integral part of jInv, due to dependency issues.
	## TODO: define the problem and fix it in a more general way for all solver with external data structures.
	## DO NOT DELETE THE CODE BELOW:
	###
    # if isa(pMis[1].pFor.Ainv,MUMPSsolver) && (pMis[1].pFor.Ainv.Ainv.ptr !=-1)
        # for ii=1:numFor
            # if any(ii.==indFors)
                # procMap[ii] = pMis[ii].pFor.Ainv.Ainv.worker
            # end
        # end
        # nextidx = p -> ( ind = find(procMap.==p);
            # if !isempty(ind);
                # ind = ind[1];
                # procMap[ind] = -1;
            # end;
            # return ind)
    # else
		nextidx = p -> (idx=i; i+=1; idx)
    # end

    matvc = zeros(length(x))
    updateRes(x) = (matvc+=x)

    @sync begin
        for p = workers()
            @async begin
                while true
                    idx = nextidx(p)
                    if isempty(idx) || idx > numFor
                        break
                    end
                    if any(idx.==indFors)
                        zi = remotecall_fetch(HessMatVec,p,x,pMis[idx],sigma,d2F[idx])
                        updateRes(zi)
                    end
                end
            end
        end
    end
    return matvc
end
