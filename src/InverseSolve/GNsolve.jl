export projGNexplicit

"""
function projGNexplicit(gc,pMis,pInv,sig,dsig,d2F,d2R,Active)

explicitly builds and solves projected normal equation.

Inputs:

	gc        - gradient
	pMis      - misfit params
	pInv      - inverse param
	sig,dsig  - current model and derivative
	d2F,d2R   - Hessians os misfit and regularizer
	Active    - indicator of active set, default false

Outputs:

	dm        - search direction
	times     - array containing time to build and solve Hessian
"""
function projGNexplicit(gc::Vector,pMis,pInv::InverseParam,sig::Vector,dsig,d2F,d2R,Active=falses(length(gc)))
	if all(Active)
		return 0*gc,[0.0;0.0]
	end
	# get Hessian of misfit
	timeBuild = @elapsed begin
	Hm = getHessian(sig,pMis,d2F)
	end;

	# build overall Hessian
	H  = dsig'*Hm*dsig + d2R

	# remove Active constraints
	timeSolve = @elapsed begin
	Hr = H[.!(Active),.!(Active)]
	gr = gc[.!(Active)]
	dm = 0*gc
	dr = -(Hr\gr)
	if any(isinf.(dr)) || any(isnan.(dr))
		dr = -(pinv(Hr)*gr)
	end
	dm[.!(Active)] = dr
	end

	# solve and return
	return dm,[timeBuild;timeSolve]
end
