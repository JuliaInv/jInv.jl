export IterativeSolver,getIterativeSolver,solveLinearSystem!

""" 
mutable struct jInv.LinearSolvers.IterativeSolver <: AbstractSolver
	
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
	isTranspose - if true, transpose(A) is provided to solver, else A is provided to solver
		      default=false, use isTranspose=true for efficiency with caution
		      note that A_mul_B! is slower than Ac_mul_B for SparseMatrixCSC

Example:

	getIterativeSolver(cg)

"""
mutable struct IterativeSolver<: AbstractSolver
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
function jInv.LinearSolvers.getIterativeSolver
	
constructs IterativeSolver

Required Input:

	IterMethod::Function   - function handle for linear solvers 
		Inputs are: (A,b,M), A is matrix, b is right hand side, M is preconditioner
			Examples: IterMethod = KrylovMethods.cg   #KrylovMethods.cg already has required API
				  IterMethod(A,b;M=M,tol=1e-1,maxIter=10,out=-1) =
				                  bicgstb(A,b,M1=M,tol=tol,maxIter=maxIter,out=out)
				  IterMethod(A,b;M=M,tol=1e-1,maxIter=10,out=-1)  = 
				                  gmres(A,b,5,M1=M,tol=tol,maxIter=maxIter,out=out)
			The keyword arguments of IterMethod for bicgstb and gmres
			will be initialized with the fields in the IterativeSolver type.
		Outputs are: (x,flag,err,iter), x is approximate solution

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
function getIterativeSolver(IterMethod::Function;PC=:ssor,maxIter=500,tol=1e-5,
					Ainv=identity,out=-1,doClear::Bool=true,nthreads::Int=4,sym=0,isTranspose=false)
 	return IterativeSolver(IterMethod,PC,maxIter,tol,Ainv,out,doClear,nthreads,sym,isTranspose,0,0,.0,.0,.0)
end

function solveLinearSystem!(A,B,X,param::IterativeSolver,doTranspose=0)
	if param.doClear
		# clear preconditioner
		clear!(param)
		param.doClear=false
	end
	
	# build preconditioner
	if param.Ainv == []
		if param.PC==:ssor
			OmInvD = 1 ./Vector(diag(A));
			x      = zeros(eltype(A),size(B,1))
			M = r -> (x[:].=0.0; param.timePC += @elapsed x=ssorPrecTrans!(A,x,r,OmInvD); return x);
			param.Ainv= M
		elseif param.PC==:jac
			OmInvD = 1 ./Vector(diag(A))
			M = r -> (param.timePC += @elapsed x=r.*OmInvD; return x); 
			param.Ainv= M
		else 
			error("Iterativesolver: preconditioner $(param.PC) not implemented.")
		end
		param.nBuildPC+=1
	end
	
	# solve systems
	y     = zeros(eltype(A),size(X,1))
	doTranspose = (param.isTranspose) ? mod(doTranspose+1,2) : doTranspose
	if hasParSpMatVec
		if (param.sym==1) ||  ((param.sym != 1) && (doTranspose == 1)) 
			Af = x -> (y[:].=0.0; param.timeMV+=@elapsed ParSpMatVec.Ac_mul_B!(one(eltype(A)),A,x,zero(eltype(A)),y,param.nthreads);  return y)
		elseif (param.sym != 1) && (doTranspose == 0)
			Af = x -> (y[:].=0.0; param.timeMV+=@elapsed ParSpMatVec.A_mul_B!(one(eltype(A)),A,x,zero(eltype(A)),y,param.nthreads);  return y)
		end
			
	else
		if (param.sym==1) ||  ((param.sym != 1) && (doTranspose == 1)) 
			Af = x -> (y[:].=0.0; param.timeMV+=@elapsed mul!(y,adjoint(A),x,one(eltype(A)),zero(eltype(A))); return y)
		elseif (param.sym != 1) && (doTranspose == 0)
			Af = x -> (y[:].=0.0; param.timeMV+=@elapsed mul!(y,A,x,one(eltype(A)),zero(eltype(A)));  return y)
		end
	end
	
	t = time_ns();
	for i=1:size(X,2)
		bi      = Vector(B[:,i]);
		X[:,i],flag,err,iter = param.IterMethod(Af,bi,M = param.Ainv,tol=param.tol,maxIter=param.maxIter,out=param.out)
		param.nIter+=iter
	end	
	param.timeSolve+=(time_ns()-t)/1e+9;
	return X, param
end # function solveLinearSystem PCGsolver
