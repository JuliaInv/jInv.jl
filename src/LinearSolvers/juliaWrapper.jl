export JuliaSolver,getJuliaSolver,copySolver,solveLinearSystem,solveLinearSystem!,setupSolver

## This did not work in 0.7. Check this in the future.
# import Base.\
# function \{T1,T2}(A::Base.SparseArrays.UMFPACK.UmfpackLU{T1},R::SparseMatrixCSC{T2})

	# n,nrhs = size(R)
	# X = zeros(promote_type(T1,T2),n,nrhs)
	# for k=1:nrhs
		# X[:,k] = A\full(vec(R[:,k]))
	# end
	# return X
# end

"""
mutable struct jInvLinearSolvers.JuliaSolver<: AbstractSolver

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

function setupSolver(A::SparseMatrixCSC,param::JuliaSolver)
	tt = time_ns();
	if param.sym==1 && isreal(A)
		param.Ainv = cholesky(A)
	elseif param.sym==2 && isreal(A)
		param.Ainv = ldlt(A)
	else
		param.Ainv = lu(A);
	end
	param.facTime+= (tt-time_ns())/1e+9; 
	param.nFac+=1
	return param;
end



function solveLinearSystem!(A::SparseMatrixCSC,B,X,param::JuliaSolver,doTranspose=0)
	if issparse(B)
		#println("");
		#@warn("jInv: Julia solvers do not support sparse RHSs for now. Check in the future");
		if length(size(B))==1
			B = Vector(B);
		else
			B = Matrix(B);
		end
	end
	if param.doClear == 1
		clear!(param)
	end
	if param.sym==0
		if doTranspose==1 && param.isTransposed==0
			clear!(param);
			A = sparse(A');
			param.isTransposed = 1;
		end
		if doTranspose==0 && param.isTransposed==1
			clear!(param);
		end
	end
	if param.Ainv == []
		param = setupSolver(A,param);
	end
	
	tt = time_ns()
	U = param.Ainv\B;
	param.solveTime+=(tt-time_ns())/1e+9; 
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
