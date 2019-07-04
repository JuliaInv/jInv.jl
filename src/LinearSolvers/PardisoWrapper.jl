export jInvPardisoSolver, getjInvPardisoSolver

"""
type jInvPardisoSolver

Fields:

	Ainv    - holds PardisoFactorization
	doClear - flag to clear factorization
	ooc     - flag for out-of-core option
	sym     - 0=unsymmetric, 1=symm. pos def, 2=general symmetric
	nFac    - number of factorizations performed
	facTime - cumulative time for factorizations
	nSolve  - number of solves
	solveTime - cumulative time for solves

Example:

	Ainv = getjInvPardisoSolver()
"""
mutable struct jInvPardisoSolver<: AbstractDirectSolver
	Ainv
	doClear::Int
	ooc::Int
	sym::Int
	nFac::Int
	facTime::Real
	nSolve::Int
	solveTime::Real
	N::Int
	nProcs::Int
end

# See if Pardiso is installed and include code for it.


if hasPardiso

	function getjInvPardisoSolver(Ainv=[],doClear=1,ooc=0,sym=11,nProcs=1)
		return jInvPardisoSolver(Ainv,doClear,ooc,sym,0,0.0,0,0.0,1,nProcs)
	end

	function factorLinearSystem!(A,param::jInvPardisoSolver)
		if param.doClear == 1
			clear!(param)
		end
        X = Vector{eltype(A)}(size(A,1))
		param.facTime+=@elapsed begin
		param.Ainv = MKLPardisoSolver()
		set_nprocs!(param.Ainv, param.nProcs)
		param.N    = size(A,1)
		param,Apard = pardisoSetup(param,A)
		pardiso(param.Ainv,X,Apard,X)
		end
		param.nFac+=1
	end

	function solveLinearSystem!(A,B,X,param::jInvPardisoSolver,doTranspose=0)
		if param.doClear == 1
			clear!(param)
		end
		if param.Ainv==[]
			param.facTime+=@elapsed begin
			param.Ainv = MKLPardisoSolver()
			set_nprocs!(param.Ainv, param.nProcs)
			param.N    = size(A,1)
			param,Apard = pardisoSetup(param,A)
			pardiso(param.Ainv,X,Apard,full(B))
			end
			param.nFac+=1
		end

		# Since pardiso requires CSR matrices as input
		# and Julia sparse matrices are CSC,
		# we need to tell pardiso to solve transposed
		# system when we want untransposed solve and
		#vice versa
		param.solveTime+=@elapsed begin
		if (param.sym in [1, 3, 11, 13]) && (doTranspose == 0)
		  set_iparm!(param.Ainv, 12, 2)
		else
		  set_iparm!(param.Ainv, 12, 0)
		end
		set_phase!(param.Ainv,33)
		set_iparm!(param.Ainv, 27, 0) #Deactivate matrix checker
		pardiso(param.Ainv,X,A,full(B))
		end
		param.nSolve+=1

		return X, param
	end # function solveLinearSystem MKLPardisoSolver

	function pardisoSetup(param::jInvPardisoSolver,A::SparseMatrixCSC)
			set_phase!(param.Ainv,12) #Perform analysis and numerical fac.
			set_matrixtype!(param.Ainv,param.sym) #Set matrix type
			pardisoinit(param.Ainv) #Initialize iparm
			set_iparm!(param.Ainv, 2, 3) #Use parallel metis ordering
			set_iparm!(param.Ainv, 27, 1) #Activate matrix checker
			if param.ooc != 0
			  set_iparm!(param.Ainv, 60, 2)
			end

			# Pardiso only needs upper triangular part of
			# symmetric matrices. We take lower triangular
			# part due to CSC vs CSR conflict.
			if param.sym in [2, 4, -2, -4, 6]
			  Apard = tril(A)
			else
			  Apard = A
			end
			return param,Apard
	end

	function clear!(param::jInvPardisoSolver)
		if param.Ainv==[]
		        return
		else
			set_phase!(param.Ainv,-1)
			# Real or complex matrix type
			if param.sym in [1, 2, -2, 11]
				A = sparse(1.0I,2,2);
				x = ones(2)
				b = ones(2)
			else
				A = sparse(one(ComplexF64)*I, 2,2)
				x = ones(ComplexF64, 2)
				b = ones(ComplexF64, 2)
			end
			pardiso(param.Ainv,x,A,b)
			param.Ainv = []
		end
	end

	function copySolver(Ainv::jInvPardisoSolver)
		return jInvPardisoSolver(Ainv.Ainv,Ainv.doClear,
		                         Ainv.ooc,Ainv.sym,0,0.0,0
		                         ,0.0,1,Ainv.nProcs);
	end

end
