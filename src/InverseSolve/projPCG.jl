export projPCG

"""
	dm = projPCG(H,g,Active,Precond,cgTol,maxIter)

	Projected Preconditioned Conjugate Gradient method for solving

		H*dm = g    subject to    dm[!Active] == 0

	Input:

		H::Function       - computes action of Hessian
		g::Vector         - right hand side
		Active            - describes active cells
		Precond::Function - preconditioner
		cgTol             - tolerance
		maxIter           - maximum number of iterations
    out               - verbosity level
"""
function projPCG(H::Function,gc::Vector,Active::BitArray,Precond::Function,cgTol::Real,maxIter::Int;out::Int=0)
#  PCG over active cells
    his        = zeros(maxIter,3)
	if norm(gc)==0
    if out >= 2
      println("projPCG finished without iterating: trivial right hand side and solution")
    end
    return zeros(eltype(gc),length(gc)),his[1,:]
  end

	delm       = zeros(eltype(gc),size(gc))
	cgiter     = 0
	resid      = - .!Active.*gc
	normResid0 = norm(resid)
	rdlast = 0.0;
	pc = 0.0;
	while true
		cgiter += 1

		his[cgiter,2]= @elapsed dc = .!Active.*Precond(resid)

		rd = dot(resid,dc)

		#  Compute conjugate direction pc.
		if cgiter == 1
			pc = dc
		else
			betak = rd / rdlast
			pc = dc + betak * pc
		end

		#  Form product Hessian*pc.

		his[cgiter,3] = @elapsed Hp = H(pc)


		Hp = .!Active.*Hp
		#  Update delm and residual.
		alphak = rd / dot(pc,Hp)
		delm   = delm + alphak*pc
		resid  = resid - alphak*Hp
		rdlast = rd

	    his[cgiter,1] = norm(resid)/normResid0;
    if out >=3
      println(@sprintf("projPCG iteration %d: relres = %.3e", cgiter, his[cgiter,1]))
    end
		if (his[cgiter,1] <= cgTol) || (cgiter == maxIter )
	  	  break
		end
	end
  if out >= 2
    if his[cgiter,1] <= cgTol
      println(@sprintf("projPCG converged at iteration %d with relres = %.3e", cgiter, his[cgiter,1]))
    else
      println(@sprintf("projPCG stopped at iteration %d with relres = %.3e", cgiter, his[cgiter,1]))
    end
  end

	return delm,his[1:cgiter,:]
end

function projPCG(gc::Vector,pMis,pInv,sig::Vector,dsig,d2F,d2R,Active;out::Int=0)

		#  Set up Hessian and preconditioner.
		Hs(x) = dsig'*HessMatVec(dsig*x,pMis,sig,d2F) + d2R*x;
		#  build preconditioner
		pInv.HesPrec.setupPrec(Hs, d2R,pInv.HesPrec.param);
		PC(x) = pInv.HesPrec.applyPrec(Hs,d2R,x,pInv.HesPrec.param);
		# call projPCG and return
		return projPCG(Hs,gc,Active,PC,pInv.pcgTol,pInv.pcgMaxIter,out=out)
end
