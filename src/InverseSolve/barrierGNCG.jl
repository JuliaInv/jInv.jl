export barrierGNCG

function shrinkStep(mc,d,low,high,out::Int=2)
fact = 1.0;
if maximum(mc + d - high) > -1e-16
	facts = (high - mc)./d;
	fact = minimum(facts[facts.>0.0]);
	d .*= fact;
end
if minimum(mc + d - low) < 1e-16
	facts = (low - mc)./d;
	fact = minimum(facts[facts.>0.0]);
	d .*= fact;
end
if fact != 1.0
	d .*= 0.8;
	if out >= 2
		println("Shrinking step by a factor of ",fact*0.8);
	end
end
return d;
end


function  barrierGNCG(mc,pInv::InverseParam,pMis;rho = 10.0,epsilon = 0.1*(pInv.boundsHigh - pInv.boundsLow), indCredit=[],dumpResults::Function = dummy,out::Int=2)

	maxIter     = pInv.maxIter      #  Max. no. iterations.
	pcgMaxIter  = pInv.pcgMaxIter   #  Max cg iters.
	pcgTol      = pInv.pcgTol       #  CG stopping tolerance.
	stepTol     = pInv.minUpdate    #  Step norm stopping tol.
	maxStep     = pInv.maxStep
	low         = pInv.boundsLow
	high        = pInv.boundsHigh
	alpha       = pInv.alpha
	mref 		= pInv.mref;

	His = getGNhis(maxIter,pcgMaxIter)
	#---------------------------------------------------------------------------
	#  Initialization.
	#---------------------------------------------------------------------------

	# Active 	= convert(BitArray,zeros(Bool,length(mc)));  # Compute zero active set

	# logBarrierReg(m,mref,M) = logBarrier(m,mref,M,low,high,epsilon);
	logBarrierReg(m,mref,M) = logBarrierSquared(m,mref,M,low,high,epsilon);
	regularizer 	 		= [pInv.regularizer;logBarrierReg];
	alpha 		 	 		= [pInv.alpha;rho];
	mref 			 		= zeros(size(mref,1),size(mref,2)+1);
	mref[:,1:size(mref,2)-1] 	= pInv.mref;


	#########################################
	## evaluate function and derivatives
	sig,dsig = pInv.modelfun(mc)
	if isempty(indCredit)
		Dc,F,dF,d2F,pMis,tMis = computeMisfit(sig,pMis,true)
	else
		Dc,F,dF,d2F,pMis,tMis,indDebit = computeMisfit(sig,pMis,true,indCredit=indCredit)
	end
	dF = dsig'*dF


	# compute regularizer
	tic()
	R,dR,d2R = computeRegularizer(regularizer,mc,mref,pInv.MInv,alpha)
	tReg = toq()
	# objective function
	Jc  = F  + R
	gc  = dF + dR

	F0 = F; J0 = Jc
	############################################################################
	##  Outer iteration.                                                        #
	############################################################################
	iter = 0
	outerFlag = -1; stepNorm=0.0

	outStr = @sprintf("%4s\t%08s\t%08s\t%08s\t%08s\n",
					  	"i.LS", "F", "R","alpha[1]","Jc/J0")
	updateHis!(0,His,Jc,norm(gc),F,Dc,R,alpha[1],0,0.0,-1,tMis,tReg)

	if out>=2; print(outStr); end
	f = open("jInv.out", "w")
	write(f, outStr)
	close(f)

	while outerFlag == -1

		iter += 1
		outStr = @sprintf("%3d.0\t%3.2e\t%3.2e\t%3.2e\t%3.2e\n",
		         iter, F, R,alpha[1],Jc/J0)
		if out>=2; print(outStr); end
		f = open("jInv.out", "a")
		write(f, outStr)
		close(f)


		#  Set up Hessian and preconditioner.
		if isempty(indCredit)
			Hs = x -> dsig'*HessMatVec(dsig*x,pMis,sig,d2F) + d2R*x;
		else
			Hs = x -> dsig'*HessMatVec(dsig*x,pMis,sig,d2F,indDebit) + d2R*x;
		end


		pInv.HesPrec.setupPrec(Hs, d2R,pInv.HesPrec.param);

		PC(x) = pInv.HesPrec.applyPrec(Hs,d2R,x,pInv.HesPrec.param);



		##  inner CG iterations.
		tic()
		delm,hisLinSol = KrylovMethods.cg(Hs,-gc, tol=pcgTol, maxIter=pcgMaxIter, M=PC,out=0)[1:2];
		His.timeLinSol[iter+1]+= toq()
		push!(His.hisLinSol,hisLinSol)


		# scale step
		if maximum(abs.(delm)) > maxStep; delm = delm./maximum(abs.(delm))*maxStep; end

		delm = shrinkStep(mc,delm,low,high,out);

		## Begin Armijo line search
		muLS = 1.0; lsIter = 1; mt = zeros(size(mc)); Jt = Jc
		while true
			mt = mc + muLS*delm

			## evaluate function
			sigt, = pInv.modelfun(mt)
			if isempty(indCredit)
				Dc,F,dF,d2F,pMis,tMis = computeMisfit(sigt,pMis,false)
			else
				Dc,F,dF,d2F,pMis,tMis,indDebit = computeMisfit(sigt,false,indCredit=indCredit)
			end
			His.timeMisfit[iter+1,:]+=tMis

			tic()
			R,dR,d2R = computeRegularizer(regularizer,mt,mref,pInv.MInv,alpha)
			His.timeReg[iter+1] += toq()
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

		#  Check stopping criteria for outer iteration.
		updateHis!(iter,His,Jc,-1.,F,Dc,R,alpha[1],0,stepNorm,lsIter,tMis,tReg)


		dumpResults(mc,Dc,iter,pInv,pMis);
		if stepNorm < stepTol
			outerFlag = 1
			break
		elseif iter >= maxIter
			break
		end

		# Evaluate gradient
		tic()
		if isempty(indCredit)
			dF = computeGradMisfit(sig,Dc,pMis)
		else
			dF = computeGradMisfit(sig,Dcp,pMis,indDebit)
		end
		His.timeGradMisfit[iter+1]+=toq()

		dF = dsig'*dF
		gc = dF + dR


		His.dJ[iter+1] = norm(gc)
	end # while outer_flag == 0

	if out>=1
		if outerFlag==-1
			println("barrierGNCG iterated maxIter=$maxIter times but reached only stepNorm of $(stepNorm) instead $(stepTol)." )
		elseif outerFlag==-2
			println("barrierGNCG stopped at iteration $iter with a line search fail.")
		elseif outerFlag==1
			println("barrierGNCG reached desired accuracy at iteration $iter.")
		end
	end


	return mc,Dc,outerFlag,His
end  # Optimization code
