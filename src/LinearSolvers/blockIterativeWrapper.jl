export BlockIterativeSolver,getBlockIterativeSolver,solveLinearSystem!

""" 
mutable struct BlockIterativeSolver
	
Fields:

	IterMethod - iterative method to apply
	PC      - symbol, (:ssor, :jac,...)
	maxIter - maximum number of iterations
	tol     - tolerance
	Ainv    - preconditioner
	out     - flag for output
	doClear - flag for deleting preconditioner
	nthreads - number of threads for spmatvecs
	sym      - 0=unsymmetric, 1=symm. pos def, 2=general symmetric
	isTranspose - if true, transpose(A) is provided to solver, else A is proved to solver
		      default=false, use isTranspose=true for efficiency with caution
		      note that A_mul_B! is slower than Ac_mul_B for SparseMatrixCSC
	

Example
getBlockIterativeSolver()
"""
mutable struct BlockIterativeSolver<: AbstractSolver
	IterMethod::Function
	PC::Symbol
	maxIter::Int
	tol::Real
	Ainv
	out::Int
	doClear::Bool
	nthreads::Int
	sym::Int
	isTranspose::Bool
	nIter::Int
	nBuildPC::Int
	timePC::Real
	timeSolve::Real
	timeMV::Real
end

"""
function jInv.LinearSolvers.getBlockIterativeSolver
	
constructs BlockIterativeSolver

Required Input:

	IterMethod::Function   - function handle for linear solvers 
		Inputs are: (A,B,M), A is matrix, B are right hand sides, M is preconditioner
			Examples: 
				  IterMethod = blockCG
				  IterMethod(A,B;M=M,X=X,tol=1e-1,maxIter=10,out=-1) =
				                  blockBiCGSTB(A,b,M1=M,X=X,tol=tol,maxIter=maxIter,out=out)
			The keyword arguments of IterMethod for blockBiCGSTB
			will be initialized with the fields in the IterativeSolver type.
		Outputs are: (X,flag,err,iter), X are approximate solutions

Optional Inputs:

	PC::Symbol     - specifies preconditioner, default:ssor
	maxIter        - maximum number of iterations, default:500
	tol            - tolerance on relative residual, default=1e-5
	Ainv           - preconditioner, default=identity
	out            - flag for output, default=-1 (no output)
	doClear        - flag for clearing the preconditioner, default=true
	nthreads       - number of threads to use for matvecs (requires ParSpMatVec.jl), default=4
	sym            - 0=unsymmetric, 1=symm. pos def, 2=general symmetric
    isTranspose    - if true, transpose(A) is provided to solver, else A is proved to solver
		              default=false, use isTranspose=true for efficiency with caution
		              note that A_mul_B! is slower than Ac_mul_B for SparseMatrixCSC
"""
function getBlockIterativeSolver(IterMethod;PC=:ssor,maxIter=10,tol=1e-4,Ainv=identity,
	out=-1,doClear::Bool=true,ortho::Bool=false,nthreads::Int=4,sym=0,isTranspose=false)
	return BlockIterativeSolver(IterMethod,PC,maxIter,tol,Ainv,out,doClear,nthreads,sym,isTranspose,0,0,.0,.0,.0)
end


function solveLinearSystem!(A,B,X,param::BlockIterativeSolver,doTranspose=0)
	if param.doClear
		# clear preconditioner
		clear!(param)
		param.doClear=false
	end
	
	
	n = size(B,1)
	nrhs = size(B,2)
	# build preconditioner The preconditioners here are symmetric anyway.
	if param.Ainv == []
		if param.PC==:ssor
			OmInvD = 1 ./Vector(diag(A));
			Xt      = zeros(n,nrhs)
			M = R -> (Xt[:].=0.0; param.timePC+=@elapsed Xt=ssorPrecTrans!(A,Xt,R,OmInvD); return Xt);
			param.Ainv= M
		elseif param.PC==:jac
			OmInvD = 1 ./Vector(diag(A))
			M = R -> (param.timePC+=@elapsed Xt=R.*OmInvD; return Xt); 
			param.Ainv= M
		else 
			error("PCGsolver: preconditioner $(param.PC) not implemented.")
		end
		param.nBuildPC+=1
	end
	
	if issparse(B)
		B = Matrix(B);
	end
	# solve systems
	Y    = zeros(n,nrhs)
	doTranspose = (param.isTranspose) ? mod(doTranspose+1,2) : doTranspose
	if hasParSpMatVec
		if (param.sym==1) ||  ((param.sym != 1) && (doTranspose == 1)) 
			Af = X -> (param.timeMV+=@elapsed ParSpMatVec.Ac_mul_B!(one(eltype(A)),A,X,zero(eltype(A)),Y,param.nthreads); return Y)
		elseif (param.sym != 1) && (doTranspose == 0)
			Af = X -> (param.timeMV+=@elapsed ParSpMatVec.A_mul_B!(one(eltype(A)),A,X,zero(eltype(A)),Y,param.nthreads); return Y)
		end
			
	else
		if (param.sym==1) ||  ((param.sym != 1) && (doTranspose == 1)) 
			Af = X -> (param.timeMV+=@elapsed mul!(Y,adjoint(A),X,one(eltype(A)),zero(eltype(A)));  return Y)
			# Ac_mul_B!(one(eltype(A)),A,X,zero(eltype(A)),Y); return Y)
		elseif (param.sym != 1) && (doTranspose == 0)
			Af = X -> (param.timeMV+=@elapsed mul!(Y,A,X,one(eltype(A)),zero(eltype(A)));  return Y)
		end
	end
	X[:].=0.0	
	t = time_ns();
	
	X,flag,err,iter = param.IterMethod(Af,B,X=X,M=param.Ainv,tol=param.tol,
										maxIter=param.maxIter,out=param.out)
	param.nIter+=iter*nrhs
	param.timeSolve+=(t-time_ns())/1e+9;
	return X, param
end # function solveLinearSystem 



