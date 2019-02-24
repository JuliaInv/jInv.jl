export GNhis, getGNhis, updateHis!

"""
mutable struct jInv.InverseSolve.GNhis
	
stored iteration history of Gauss-Newton methods and variants

Fields:

Jc             - objective function values
dJ             - norm of (projected) gradients
F              - misfits
Dc             - data
Rc             - value of regularizer
alphas         - regularization parameters
Active         - active sets
stepNorm       - norm of step
lsIter         - nuber of line search iterations
timeMisfit     - time to evaluate misfits
timeReg        - time to evaluate regularizer
hisLinSol      - history of linear solver (if applicable)
timeLinSol     - time for linear solves
timeGradMisfit - time for gradient computation

"""
mutable struct GNhis
	Jc::Array
	dJ::Array
	F::Array
	Dc::Array
	Rc::Array
	alphas::Array
	Active::Array
	stepNorm::Array
	lsIter::Array
	timeMisfit::Array
	timeReg::Array
	hisLinSol::Array
	timeLinSol::Array
	timeGradMisfit::Array
end

"""
function jInv.InverseSolve.getGNhis(maxIter,maxIterCG)
	
constructs GNhis

Input:

maxIter
maxIterCG

"""
function getGNhis(maxIter,maxIterCG)
	Jc = zeros(maxIter+1)
	dJ = zeros(maxIter+1)
	F  = zeros(maxIter+1)
	Dc = []
	Rc = zeros(maxIter+1)
	alphas = zeros(maxIter+1)
	Active = zeros(maxIter+1)
	stepNorm = zeros(maxIter+1)
	lsIter = zeros(Int,maxIter+1)
	timeMisfit = zeros(maxIter+1,4)
	timeReg = zeros(maxIter+1)
	hisLinSol = []
	timeLinSol = zeros(maxIter+1,1)
	timeGradMisfit = zeros(maxIter+1,2)

	return GNhis(Jc,dJ,F,Dc,Rc,alphas,Active,stepNorm,lsIter,timeMisfit,timeReg,hisLinSol,timeLinSol,timeGradMisfit)
end

function updateHis!(iter::Int64,His::GNhis,Jc::Real,dJ::Real,Fc,Dc,Rc::Real,alpha::Real,nActive::Int64,stepNorm::Real,lsIter::Int,timeMisfit::Vector,timeReg::Real)
	His.Jc[iter+1]            = Jc
	His.dJ[iter+1]            = dJ
	His.F[iter+1]             = Fc
	push!(His.Dc,Dc)
	His.Rc[iter+1]            = Rc
	His.alphas[iter+1]        = alpha
	His.Active[iter+1]        = nActive
	His.stepNorm[iter+1]      = stepNorm
	His.lsIter[iter+1]        = lsIter
	His.timeMisfit[iter+1,:] += timeMisfit
	His.timeReg[iter+1]      += timeReg[]
end