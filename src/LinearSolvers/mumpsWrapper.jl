export MUMPSsolver, getMUMPSsolver,copySolver

"""
type MUMPSsolver

Fields:

	Ainv    - holds MUMPSfactorization
	doClear - flag to clear factorization
	ooc     - flag for out-of-core option
	sym     - 0=unsymmetric, 1=symm. pos def, 2=general symmetric
	nFac    - number of factorizations performed
	facTime - cumulative time for factorizations
	nSolve  - number of solves
	solveTime - cumnulative time for solves

Example:

	Ainv = getMUMPSsolver()
"""
type MUMPSsolver<: AbstractDirectSolver
	Ainv
	doClear::Int
	ooc::Int
	sym::Int
	nFac::Int
	facTime::Real
	nSolve::Int
	solveTime::Real
end

# See if MUMPS solver is installed and include code for it.


if hasMUMPS

	function getMUMPSsolver(Ainv=[],doClear=1,ooc=0,sym=0)
		return MUMPSsolver(Ainv,doClear,ooc,sym,0,0.0,0,0.)
	end

    function factorLinearSystem!(A,param::MUMPSsolver)
		if param.doClear == 1
			clear!(param)
		end

		tic()
		param.Ainv = factorMUMPS(A, param.sym, param.ooc)
		param.facTime+=toq()
		param.nFac+=1
	end

	function solveLinearSystem!(A,B,X,param::MUMPSsolver,doTranspose::Int=0)
		if param.doClear == 1
			clear!(param)
		end

		if param.Ainv==[]
			tic()
			param.Ainv = factorMUMPS(A, param.sym, param.ooc)
			param.facTime+=toq()
			param.nFac+=1
		end

		tic()

		applyMUMPS!(param.Ainv, B,X,doTranspose)

		param.solveTime+=toq()
		param.nSolve+=1

		return X, param
	end # function solveLinearSystem MUMPSsolver

	function clear!(M::MUMPSsolver)
		if isa(M.Ainv,MUMPSfactorization)
			destroyMUMPS(M.Ainv)
			M.Ainv = []
		end
	end

	function copySolver(Ainv::MUMPSsolver)
		return MUMPSsolver([],Ainv.doClear,Ainv.ooc,Ainv.sym,0,0.0,0,0.);
	end

	hasMUMPS=true
end
