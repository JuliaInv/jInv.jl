export HessianPreconditioner,getSSORCGRegularizationPreconditioner,getSSORRegularizationPreconditioner,getExactSolveRegularizationPreconditioner
export getEmptyRegularizationPreconditioner
############ General Hessian Preconditioner:
mutable struct HessianPreconditioner
	param
	applyPrec::Function # a function with parameters: Hs::function, d2R::SparseMatrixCSC,v::Vector,param::Any
	setupPrec::Function # a function with parameters: Hs::function, d2R::SparseMatrixCSC,v::Vector,param::Any
end
######################################################################################################################
############ SSORCG: A preconditioner that inverts the regularization matrix by SSOR - CG ############################
###################################################################################################################### 

mutable struct SSORCGParam
	diagonal	::Vector
	auxVec		::Vector
	omega		::Float64
	tol			::Float64
	maxCGIter	::Int64
end

function getSSORCGRegularizationPreconditioner(omega::Float64=1.0,tol::Float64=1e-2,maxCGIter::Int64=100)
	return HessianPreconditioner(SSORCGParam([],[],omega,tol,maxCGIter),applySSORCG,setupSSORCG);
end

function applySSORCG(Hs::Function, d2R::SparseMatrixCSC,v::Vector,param)
	aux = param.auxVec;
	SSOR(r) = (aux[:].=0.0; xt=ssorPrecTrans!(d2R,aux,r,param.diagonal); return aux);
	x = KrylovMethods.cg(d2R,v,tol=param.tol,maxIter=param.maxCGIter,M=SSOR,out=-1)[1]
	return x;
end

function setupSSORCG(Hs::Function, d2R::SparseMatrixCSC,param)
	param.diagonal = param.omega./diag(d2R);
	param.auxVec   = zeros(size(d2R,2));
	return;
end

######################################################################################################################
############ SSOR: A preconditioner that inverts the regularization matrix by SSOR only ##############################
###################################################################################################################### 

mutable struct SSORParam
	omega		::Float64
	tol			::Float64
	maxIter	::Int64
end

function getSSORRegularizationPreconditioner(omega::Float64=1.0,tol::Float64=1e-2,maxIter::Int64=100)
	return HessianPreconditioner(SSORParam(omega,tol,maxIter),applySSOR,setupSSOR);
end

function applySSOR(Hs::Function, d2R::SparseMatrixCSC,v::Vector,param)
	x = KrylovMethods.ssor(d2R,copy(v);x=[],tol=param.tol,maxIter=param.maxIter,omega=param.omega,out=-1,storeInterm=false)[1];
	return x;
end

function setupSSOR(Hs::Function, d2R::SparseMatrixCSC,param)
	return;
end

######################################################################################################################
############ ExactSolveRegularization: A preconditioner that inverts the regularization matrix by a backslash operator
######################################################################################################################

function getExactSolveRegularizationPreconditioner()
	return HessianPreconditioner([],applyExactRegSolve,setupExactRegSolve);
end

function applyExactRegSolve(Hs::Function, d2R::SparseMatrixCSC,v::Vector,param)
	return d2R\v;
end

function setupExactRegSolve(Hs::Function, d2R::SparseMatrixCSC,param)
	return;
end

########################################################################################
function getEmptyRegularizationPreconditioner()
	return HessianPreconditioner([],applyEmptyRegSolve,setupEmptyRegSolve);
end

function applyEmptyRegSolve(Hs::Function, d2R::SparseMatrixCSC,v::Vector,param)
	return copy(v);
end

function setupEmptyRegSolve(Hs::Function, d2R::SparseMatrixCSC,param)
	return;
end





