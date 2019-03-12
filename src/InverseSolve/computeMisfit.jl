export computeMisfit


"""
function jInv.InverseSolve.computeMisfit(...)

Computes misfit for PDE parameter estimation problem. There are several
ways to use computeMisfit, some of which are parallelized.

Inputs:

    sig  - current model
    pMis - description of misfit term (MisfitParam, Array{MisfitParam}, or Array{Future})

Optional arguments:

    doDerivative - flag for computing derivatives (default: true)
    doClear      - flag to clear memory after computation (default: false)

Output:

   Dc    - current data
   F     - current misfit
   dF    - gradient
   d2F   - Hessian of misfit
   pMis  - modified misfit param (e.g., with PDE factorizations)
   times - some runtime statistics
"""
function computeMisfit(sig,
                       pMis::MisfitParam,doDerivative::Bool=true, doClear::Bool=false;
                       printProgress=false)
    if printProgress
        error("Print progress only works with multiple pFors")
    end
#=
 computeMisfit for a single forward problem. Everything is stored in memory on the node executing this function.
=#

    times = zeros(4)
    sigma,dsigma = pMis.modelfun(sig)
    times[1] = @elapsed   sigmaloc = interpGlobalToLocal(sigma,pMis.gloc.PForInv,pMis.gloc.sigmaBackground);    
    times[2] = @elapsed   Dc,pMis.pFor  = getData(sigmaloc,pMis.pFor)      # fwd model to get predicted data
    times[3] = @elapsed   F,dF,d2F = pMis.misfit(Dc,pMis.dobs,pMis.Wd)
    if doDerivative
        times[4] = @elapsed dF = dsigma'*interpLocalToGlobal(getSensTMatVec(dF,sigmaloc,pMis.pFor),pMis.gloc.PForInv)
    end

    if doClear; clear!(pMis.pFor.Ainv); end
    return Dc,F,dF,d2F,pMis,times
end


function computeMisfit(sigmaRef::RemoteChannel,
                        pMisRef::RemoteChannel,
				      dFRef::RemoteChannel,
                  doDerivative,doClear::Bool=false)
#=
 computeMisfit for single forward problem

 Note: model (including interpolation matrix) and forward problems are RemoteRefs
=#

    rrlocs = [ pMisRef.where  dFRef.where]
    if !all(rrlocs .== myid())
        warn("computeMisfit: Problem on worker $(myid()) not all remote refs are stored here, but rrlocs=$rrlocs")
    end

    sigma = fetch(sigmaRef)
    pMis  = take!(pMisRef)

    Dc,F,dFi,d2F,pMis,times = computeMisfit(sigma,pMis,doDerivative,doClear)

    put!(pMisRef,pMis)
    # add to gradient
    if doDerivative
        dF = take!(dFRef)
        put!(dFRef,dF += dFi)
    end
    # put predicted data and d2F into remote refs (no need to communicate them)
    Dc  = remotecall(identity,myid(),Dc)
    d2F = remotecall(identity,myid(),d2F)

    return Dc,F,d2F,times
end


function computeMisfit(sigma,
	pMisRefs::Array{RemoteChannel,1},
	doDerivative::Bool=true;
	indCredit::AbstractVector=1:length(pMisRefs),
    printProgress::Bool=false)
#=
computeMisfit for multiple forward problems

This method runs in parallel (iff nworkers()> 1 )

Note: ForwardProblems and Mesh-2-Mesh Interpolation are RemoteRefs
    (i.e. they are stored in memory of a particular worker).
=#

    n = 1

	F   = 0.0
	dF  = (doDerivative) ? zeros(length(sigma)) : []
	d2F = Array{Any}(undef, length(pMisRefs));
	Dc  = Array{Future}(undef,size(pMisRefs))

	indDebit = []
	updateRes(Fi,idx) = (F+=Fi;push!(indDebit,idx))
	updateDF(x) = (dF+=x)

    workerList = []
    for k=indCredit
        push!(workerList,pMisRefs[k].where)
    end
    workerList = unique(workerList)
    sigRef = Array{RemoteChannel}(undef,maximum(workers()))
	dFiRef = Array{RemoteChannel}(undef,maximum(workers()))

	times = zeros(4);
	updateTimes(tt) = (times+=tt)

	@sync begin
		for p=workerList
			@async begin
				# communicate model and allocate RemoteRef for gradient
				sigRef[p] = initRemoteChannel(identity,p,sigma)   # send conductivity to workers
				dFiRef[p] = initRemoteChannel(zeros,p,length(sigma)) # get remote Ref to part of gradient
				# solve forward problems
				for idx in indCredit
					if pMisRefs[idx].where==p
						Dc[idx],Fi,d2F[idx],tt = remotecall_fetch(computeMisfit,p,sigRef[p],pMisRefs[idx],dFiRef[p],doDerivative)
						updateRes(Fi,idx)
						updateTimes(tt)
                        if printProgress && ((length(indDebit)/length(indCredit)) > n*0.1)
                            if doDerivative
                                println("Misfit and gradients computed for $(10*n)% of forward problems")
                            else
                                println("Misfit and gradients computed for $(10*n)% of forward problems")
                            end
                            n += 1
                        end
					end
				end

				# sum up gradients
				if doDerivative
					updateDF(fetch(dFiRef[p]))
				end
			end
		end
	end
	return Dc,F,dF,d2F,pMisRefs,times,indDebit
end


function computeMisfit(sigma,pMis::Array,doDerivative::Bool=true,indCredit=collect(1:length(pMis));
                       printProgress=false)
	#
	#	computeMisfit for multiple forward problems
	#
	#	This method runs in parallel (iff nworkers()> 1 )
	#
	#	Note: ForwardProblems and Mesh-2-Mesh Interpolation are stored on the main processor
	#		  and then sent to a particular worker, which returns an updated pFor.
	#
	numFor   = length(pMis)
 	F        = 0.0
    dF       = (doDerivative) ? zeros(length(sigma)) : []
 	d2F      = Array{Any}(undef,numFor)
 	Dc       = Array{Any}(undef,numFor)
	indDebit = []

	# draw next problem to be solved
	nextidx() = (idx = (isempty(indCredit)) ? -1 : pop!(indCredit))

 	updateRes(Fi,dFi,idx) = (F+=Fi; dF= (doDerivative) ? dF+dFi : []; push!(indDebit,idx))

	times = zeros(4);
	updateTimes(tt) = (times+=tt)

 	@sync begin
 		for p = workers()
 				@async begin
 					while true
 						idx = nextidx()
 						if idx == -1
 							break
 						end
 							Dc[idx],Fi,dFi,d2F[idx],pMis[idx],tt = remotecall_fetch(computeMisfit,p,sigma,pMis[idx],doDerivative)
 							updateRes(Fi,dFi,idx)
							updateTimes(tt)
 					end
 				end
 		end
 	end

 	return Dc,F,dF,d2F,pMis,times,indDebit
 end
