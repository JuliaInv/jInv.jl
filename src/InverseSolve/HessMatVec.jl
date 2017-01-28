# This file contains implementations of HessMatVec for different types

export HessMatVec

# Hessian is diagonal matrix and stored as a vector.
function HessMatVec(d2F::Array{Float64,1}, x::Array{Float64,1})
	return d2F .* x
end

function HessMatVec(d2F::Array{Complex128,1}, x::Array{Complex128,1})
	return complex(real(d2F) .* real(x), imag(d2F) .* imag(x))
end

# Hessian is stored as a sparse matrix.
function HessMatVec(d2F::SparseMatrixCSC{Float64}, x::Array{Float64,1})
	return d2F * x
end

function HessMatVec(d2F::SparseMatrixCSC{Complex128}, x::Array{Complex128,1})
	return complex(real(d2F) * real(x), imag(d2F) * imag(x))
end

function HessMatVec(d2F::SparseMatrixCSC{Float64}, x::Array{Complex128,1})
   return real2complex( d2F * complex2real(x))
end



function HessMatVec(x,
                    pMis::MisfitParam,
                    sig,  # conductivity on inv mesh (active cells only)
                    d2F)
	#
	#	JTx = HessMatVec(x,pFor,mc,model,d2F)
	#
	#	Hessian * Vector for one forward problem
	#
	#	Note: all variables have to be in memory of the worker executing this method
	#
	try
		# compute local model (for multiphysics)
		sigma,dsigma = pMis.modelfun(sig)

		sigmaloc = interpGlobalToLocal(sigma,pMis.gloc.PForInv,pMis.gloc.sigmaBackground)
		xloc     = interpGlobalToLocal(dsigma*x,pMis.gloc.PForInv)
		Jx       = vec(getSensMatVec(xloc,sigmaloc,pMis.pFor))   # Jx = J*(dsig/dm)*x
		Jx       = HessMatVec(d2F,Jx) # chain rule for misfit function
		JTxloc   = getSensTMatVec(Jx,sigmaloc,pMis.pFor)
		JTx      = dsigma'*interpLocalToGlobal(JTxloc,pMis.gloc.PForInv) # = (dsig/dm)'*J'*d2F*J*(dsig/dm)*x
		return JTx
	catch err
		if isa(err,InterruptException)
			return -1
		else
			throw(err)
		end
	end
end # function HessMatVec


function HessMatVec(xRef::Future,
                    pMisRef::RemoteChannel,
                    sigmaRef::Future,
                    d2FRef::Future,
					mvRef::RemoteChannel)
	#
	#	JTx = HessMatVec(x,pFor,mc,model,d2F)
	#
	#	Hessian * Vector for one forward problem
	#
	#	Note: model and forward problem are represented as RemoteReferences.
	#
	#	!! make sure that everything is stored on this worker to avoid communication !!
	#

	rrlocs = [xRef.where pMisRef.where sigmaRef.where d2FRef.where mvRef.where]
	if !all(rrlocs .== myid())
		warn("HessMatVec: Problem on worker $(myid()) not all remote refs are stored here, but rrlocs=$rrlocs")
	end

	# fetching and taking: should be no-ops as all RemoteRefs live on my worker
	tic()
	x     = fetch(xRef)
	sigma = fetch(sigmaRef)
	pMis  = take!(pMisRef)
	d2F   = fetch(d2FRef)
	commTime = toq()

	# compute HessMatVec
	tic()
	mvi  = HessMatVec(x,pMis,sigma,d2F)
	compTime = toq()

	# putting: should be no ops.
	tic()
	mv   = take!(mvRef)
	put!(mvRef,mv+mvi)
	put!(pMisRef,pMis)
	commTime +=toq()

#	if commTime/compTime > 1.0
#		warn("HessMatVec:  Communication time larger than computation time! commTime/compTime = $(commTime/compTime)")
#	end
	return true,commTime,compTime
end


function HessMatVec(x,
					pMisRefs::Array{RemoteChannel,1},
					sigma,
					d2F,
					indFors=1:length(pMisRefs))
  #
	#	JTx = HessMatVec(x,pFor,mc,model,d2F)
	#
	#	Hessian * Vector for multiple forward problem
	#
	#	Note: models and forward problems are represented as RemoteReferences.
	#

	tic()

	# find out which workers are involved
	workerList = []
	for k=1:length(pMisRefs)
		push!(workerList,pMisRefs[k].where)
	end
	workerList = unique(workerList)
	# send sigma to all workers
	sigmaRef = Array(Future,maximum(workers()))
	yRef = Array(Future,maximum(workers()))
	zRef = Array(RemoteChannel,maximum(workers()))
	z = zeros(length(x))
	commTime = 0.0
	compTime = 0.0
	updateTimes(c1,c2) = (commTime+=c1; compTime+=c2)
	updateMV(x) = (z+=x)

	tic()
	c2 = toq()
	updateTimes(0.0,c2)

	tic()
	@sync begin
		for p = workerList
			@async begin
				# send model and vector to worker
				tic()
				sigmaRef[p] = remotecall(identity,p,sigma)
				yRef[p]     = remotecall(identity,p,x)
				zRef[p]     = initRemoteChannel(identity,p,zeros(length(x)))
				c1 = toq()
				updateTimes(c1,0.0)

				# do the actual computation
				for idx=indFors
					if pMisRefs[idx].where==p
						b,c1,c2 = remotecall_fetch(HessMatVec,p,yRef[p],pMisRefs[idx],sigmaRef[p],d2F[idx],zRef[p])
						updateTimes(c1,c2)
					end
				end

				# get back result (timing is a mixture of addition and computation)
				tic()
				updateMV(fetch(zRef[p]))
				c1 = toq()
				updateTimes(c1,0.0)
			end
		end
	end
	chkTime =toq()

	tic()
	matvc = z
	c2 = toq()
	updateTimes(0.0,c2)

	totalTime = toq()
#	@printf "Timings for HessMatVec - total: %1.3f, communication: %1.3f, computation: %1.3f\n" totalTime commTime compTime
	return matvc
end  # function HessMatVec


function HessMatVec(x,pMis::Array{MisfitParam},sigma,d2F,indFors=1:length(pMis))
	#
	#	JTx = HessMatVec(x,pFor,mc,model,d2F)
	#
	#	Hessian * Vector for multiple forward problem
	#
	#	Note: models, interpolations and forward problems are stored on main process and then sent to workers
	#
	#	!! this method may lead to more communication than the ones above !!
	#
	numFor = length(pMis);

	# get process map
	i      = 1; nextidx = p -> i;
	procMap = zeros(Int64,numFor)
	if isa(pMis[1].pFor.Ainv,MUMPSsolver) && (pMis[1].pFor.Ainv.Ainv.ptr !=-1)
		for ii=1:numFor
			if any(ii.==indFors)
				procMap[ii] = pMis[ii].pFor.Ainv.Ainv.worker
			end
		end
		nextidx = p -> ( ind = find(procMap.==p);
			if !isempty(ind);
				ind = ind[1];
				procMap[ind] = -1;
			end;
			return ind)
	else
		nextidx = p -> (idx=i; i+=1; idx)
	end

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
