export JuliaSolver,getJuliaSolver,copySolver
using SuiteSparse

import Base.\
function \(A::SuiteSparse.UMFPACK.UmfpackLU{T1},R::SparseMatrixCSC{T2}) where {T1,T2}

	n,nrhs = size(R)
	X = zeros(promote_type(T1,T2),n,nrhs)
	for k=1:nrhs
		X[:,k] = A\Vector(vec(R[:,k]))
	end
	return X
end

"""
mutable struct jInv.LinearSolvers.JuliaSolver<: AbstractSolver

Fields:

	Ainv         - holds factorization (LU or Cholesky)
	sym          - 0=unsymmetric, 1=symm. pos def, 2=general symmetric
	isTransposed - flag whether A comes transposed or not
	doClear      - flag to clear factorization
	facTime      - cumulative time for factorizations
	nSolve       - number of solves
	solveTime    - cumnulative time for solves
	nFac         - number of factorizations performed

Example:

	Ainv = getJuliaSolver()
"""
mutable struct JuliaSolver<: AbstractSolver
	Ainv
	sym::Int # 0 = unsymmetric, 1 = symmetric s.t A = A';
	isTransposed::Int
	doClear::Int
	facTime::Real
	nSolve::Int
	solveTime::Real
	nFac::Int
end


"""
function jInv.LinearSolvers.getJuliaSolver

Constructor for JuliaSolver

Optional Keyword Arguments

	Ainv = []
	sym = 0
	isTransposed = 0
	doClear = 0
"""
function getJuliaSolver(;Ainv = [],sym = 0,isTransposed = 0, doClear = 0)
	return JuliaSolver(Ainv,sym,isTransposed,doClear,0.0,0,0.0,0);
end

solveLinearSystem(A,B,param::JuliaSolver,doTranspose::Int=0) = solveLinearSystem!(A,B,[],param,doTranspose)

function solveLinearSystem!(A::SparseMatrixCSC,B,X,param::JuliaSolver,doTranspose=0)
	if param.doClear == 1
		clear!(param)
	end
	if param.sym==0
		if doTranspose==1 && param.isTransposed==0
			clear!(param);
			A = A';
			param.isTransposed = 1;
		end
		if doTranspose==0 && param.isTransposed==1
			clear!(param);
		end
	end
	if param.Ainv == []
		tic()
		if param.sym==1 && isreal(A)
			param.Ainv = chol(A)
		elseif param.sym==2 && isreal(A)
			param.Ainv = ldlt(A)
		else
			if param.sym!=0 && !isreal(A)
				warn("jInv.JuliaSolver: using lufact for complex matrix")
			end

		  param.Ainv = lu(A)
		end
		param.facTime+=toq()
		param.nFac+=1
	end

	tic()
	U = param.Ainv\B;
	param.solveTime+=toq()
	param.nSolve+=1

	return U, param
end # function solveLinearSystem

function clear!(param::JuliaSolver)
	param.Ainv = [];
	param.isTransposed = 0;
	param.doClear = 0;
end

function copySolver(Ainv::JuliaSolver)
	return getJuliaSolver(sym = Ainv.sym,doClear = Ainv.doClear);
end
