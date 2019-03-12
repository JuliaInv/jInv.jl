export projSD, projSDhis

mutable struct  projSDhis
	Jc::Array
	dJ::Array
	F::Array
	Dc::Array
	Rc::Array
	alphas::Array
	Active::Array
	stepNorm::Array
	lsIter::Array
	timeMisfit::Array
	timeReg::Array
	timeGradMisfit::Array
end

function getProjSDhis(maxIter)
	Jc             = zeros(maxIter+1)
	dJ             = zeros(maxIter+1)
	F              = zeros(maxIter+1)
	Dc             = []
	Rc             = zeros(maxIter+1)
	alphas         = zeros(maxIter+1)
	Active         = zeros(maxIter+1)
	stepNorm       = zeros(maxIter+1)
	lsIter         = zeros(Int,maxIter+1)
	timeMisfit     = zeros(maxIter+1,4)
	timeReg        = zeros(maxIter+1)
	timeGradMisfit = zeros(maxIter+1,2)

	return projSDhis(Jc,dJ,F,Dc,Rc,alphas,Active,stepNorm,lsIter,timeMisfit,timeReg,timeGradMisfit)
end

function updateHis!(iter::Int64,His::projSDhis,Jc::Real,dJ::Real,Fc,Dc,Rc::Real,alpha::Real,
					nActive::Int64,stepNorm::Real,lsIter::Int,timeMisfit::Vector,timeReg::Real)
	His.Jc[iter+1]            = Jc
	His.F[iter+1]             = Fc
	push!(His.Dc,Dc)
	His.Rc[iter+1]            = Rc
	His.alphas[iter+1]        = alpha
	His.Active[iter+1]        = nActive
	His.stepNorm[iter+1]      = stepNorm
	His.timeMisfit[iter+1,:] += timeMisfit
	His.timeReg[iter+1]      += timeReg
end


"""
	mc,Dc,outerFlag = projSD(mc,pInv::InverseParam,pMis, indFor = [], dumpResults::Function = dummy)

	(Projected) Steepest Descent method for solving

		min_x misfit(x) + regularizer(x) subject to  x in C

	where C is a convex set and a projection operator proj(x) needs to be provided.

	Input:

		mc::Vector          - intial guess for model
		pInv::InverseParam  - parameter for inversion
		pMis                - misfit terms
		indCredit           - indices of forward problems to work on
		dumpResults			- A function pointer for saving the results throughout the iterations.
							- We assume that dumpResults is dumpResults(mc,Dc,iter,pInv,pMis),
							- where mc is the recovered model, Dc is the predicted data.
							- If dumpResults is not given, nothing is done (dummy() is called).
		out::Int            - flag for output (-1: no output, 1: final status, 2: residual norm at each iteration)

	Output:
		mc                  - final model
		Dc                  - data
		outerFlag           - flag for convergence
		His                 - iteration history

"""
function  projSD(mc,pInv::InverseParam,pMis; proj=x->min.(max.(x,pInv.boundsLow),pInv.boundsHigh),
	indCredit=[], dumpResults::Function = dummy,out::Int=2)

	maxIter     = pInv.maxIter      #  Max. no. iterations.
	stepTol     = pInv.minUpdate    #  Step norm stopping tol.
	maxStep     = pInv.maxStep
	alpha       = pInv.alpha
	low         = pInv.boundsLow
	high        = pInv.boundsHigh

	His = getProjSDhis(maxIter)
	#---------------------------------------------------------------------------
	#  Initialization.
	#---------------------------------------------------------------------------
	mc = proj(mc)

	Active = (mc .<=low) .| (mc.>=high)  # Compute active set


	## evaluate function and derivatives
	sig,dsig = pInv.modelfun(mc)
	if isempty(indCredit)
		Dc,F,dF,d2F,pMis,tMis = computeMisfit(sig,pMis,true)
	else
		Dc,F,dF,d2F,pMis,tMis,indDebit = computeMisfit(sig,pMis,true,indCredit=indCredit)
	end
	dF = dsig'*dF


	# compute regularizer
	tReg = @elapsed R,dR,d2R = computeRegularizer(pInv.regularizer,mc,pInv.mref,pInv.MInv,alpha)

	# objective function
	Jc  = F  + R
	gc  = dF + dR

	F0 = F; J0 = Jc
	############################################################################
	##  Outer iteration.                                                        #
	############################################################################
	iter = 0
	outerFlag = -1; stepNorm=0.0

	outStr = @sprintf("%4s\t%08s\t%08s\t%08s\t%08s\t%08s\n",
					  	"i.LS", "F", "R","alpha[1]","Jc/J0","#Active")
	updateHis!(0,His,Jc,norm(projGrad(gc,mc,low,high)),F,Dc,R,alpha[1],count(!iszero, Active),0.0,-1,tMis,tReg)

	if out>=2; print(outStr); end
	f = open("projSD.out", "w")
	write(f, outStr)
	close(f)

	while outerFlag == -1

		iter += 1
		outStr = @sprintf("%3d.0\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3d\n",
		         iter, F, R,alpha[1],Jc/J0,count(!iszero, Active))
		if out>=2; print(outStr); end
		f = open("jInv.out", "a")
		write(f, outStr)
		close(f)


		# scale step
		if maximum(abs.(gc)) > maxStep; gc = gc./maximum(abs.(gc))*maxStep; end

		## Begin projected Armijo line search
		muLS = 1; lsIter = 1; mt = zeros(size(mc)); Jt = Jc
		while true
			mt = proj(mc - muLS*gc)
			## evaluate function
			sigt, = pInv.modelfun(mt)
			if isempty(indCredit)
				Dc,F,dF,d2F,pMis,tMis = computeMisfit(sigt,pMis,false)
			else
				Dc,F,dF,d2F,pMis,tMis,indDebit = computeMisfit(sigt,false,indCredit=indCredit)
			end
			His.timeMisfit[iter+1,:]+=tMis

			
			His.timeReg[iter+1] += @elapsed R,dR,d2R = computeRegularizer(pInv.regularizer,mt,pInv.mref,pInv.MInv,alpha)

			# objective function
			Jt  = F  + R
			if out>=2;
				println(@sprintf( "   .%d\t%3.2e\t%3.2e\t\t\t%3.2e",
			           lsIter, F,       R,       Jt/J0))
			end

			if Jt < Jc
			    break
			end
			muLS /=2; lsIter += 1
			if lsIter > 6
			    outerFlag = -2
				break
			end
		end
		## End Line search

		## Check for termination
		stepNorm = norm(mt-mc,Inf)
		mc = mt
		Jc = Jt

		sig, dsig = pInv.modelfun(mc)

		Active = (mc .<=low) .| (mc.>=high)  # Compute active set

		#  Check stopping criteria for outer iteration.
		updateHis!(iter,His,Jc,-1,F,Dc,R,alpha[1],count(!iszero, Active),stepNorm,lsIter,tMis,tReg)

		dumpResults(mc,Dc,iter,pInv,pMis);
		if stepNorm < stepTol
			outerFlag = 1
			break
		elseif iter >= maxIter
			break
		end
		# Evaluate gradient
		t = time_ns();
		if isempty(indCredit)
			dF = computeGradMisfit(sig,Dc,pMis)
		else
			dF = computeGradMisfit(sig,Dcp,pMis,indDebit)
		end
		His.timeGradMisfit[iter+1]+=(time_ns() - t)/1e+9;

		dF = dsig'*dF
		gc = dF + dR
		His.dJ[iter+1] = norm(projGrad(gc,mc,low,high))
	end # while outer_flag == 0

	if out>=1
		if outerFlag==-1
			println("projSD iterated maxIter=$maxIter times but reached only stepNorm of $(stepNorm) instead $(stepTol)." )
		elseif outerFlag==-2
			println("projSD stopped at iteration $iter with a line search fail.")
		elseif outerFlag==1
			println("projSD reached desired accuracy at iteration $iter.")
		end
	end

	return mc,Dc,outerFlag,His
end  # Optimization code
