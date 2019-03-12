export barrierGNCG

function  barrierGN(mc,pInv::InverseParam,pMis;rho = 10.0,epsilon = 0.1*(pInv.boundsHigh - pInv.boundsLow), indCredit=[],dumpResults::Function = dummy,out::Int=2,solveGN::Function=projPCG)
	low         = pInv.boundsLow .+ 1e-6; # this is not to get confused with the bounds of projGN.
	high        = pInv.boundsHigh .- 1e-6;
	#---------------------------------------------------------------------------
	#  Initialization.
	#---------------------------------------------------------------------------
	
	# logBarrierReg(m,mref,M) = logBarrier(m,mref,M,low,high,epsilon);
	logBarrierReg(m,mref,M)     = logBarrierSquared(m,mref,M,low,high,epsilon);
	pInv.regularizer 	 		= [pInv.regularizer;logBarrierReg];
	pInv.alpha 		 	 		= [pInv.alpha;rho];
	mref1 						= pInv.mref;
	mref2 			 			= zeros(size(mref1,1),size(mref1,2)+1);
	mref2[:,1:size(mref1,2)] 	= mref1;
	pInv.mref = mref2;
	
	mc,Dc,outerFlag,His = projGN(mc,pInv,pMis,indCredit=indCredit, dumpResults = dumpResults,out=out,solveGN=solveGN);
	
	
	if length(pInv.regularizer) > 2
		pInv.regularizer = pInv.regularizer[1:end-1];
		pInv.alpha = pInv.alpha[1:end-1];
	else # if we choose only one variable in the array, julia does not make it a scalar.
		pInv.regularizer = pInv.regularizer[1];
		pInv.alpha = pInv.alpha[1];
	end
	pInv.mref = mref1;
	
	
	return mc,Dc,outerFlag,His
end  # Optimization code
barrierGNCG = barrierGN;