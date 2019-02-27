using MAT # include only on main worker
using Test
using DivSigGrad
using jInv.InverseSolve
using jInv.Mesh
using jInv.LinearSolvers
using jInv.ForwardShare
using jInv.Utils

using EikonalInv

plotting = false;
if plotting
	using jInvVis
	using PyPlot
end

@everywhere begin
	include("DriversDC/velToConductMod.jl");
end


include("DriversEikonalInv/readModelAndGenerateMeshMref.jl");
include("DriversEikonalInv/prepareTravelTimeDataFiles.jl");
include("DriversEikonalInv/setupTravelTimeTomography.jl");
include("DriversDC/getSourcesAndReceivers.jl")
include("DriversDC/setupDCAndGetData.jl");


invertDC  = true;
invertEik = true;
invertJoint = invertDC & invertEik

	
matfile   = matread("../3Dseg12812864.mat")
mtrue     = matfile["VELc"]
n         = [64;64;32]
nfwd      = [32;32;16]
domain    = [0;13.5;0;13.5;0;4.2]
mtrue*=1e-3; ## mtrue is in km per sec.

###########################################################################################################################
################ Parameters For Inversion #################################################################################
###########################################################################################################################
(mtrue,MinvNodal,mref,~,~)  = readModelAndGenerateMeshMref(mtrue,0,n,domain);
Minv    					= getRegularMesh(MinvNodal.domain,collect(size(mtrue)));
HesPrec     = getSSORRegularizationPreconditioner(1.0,1e-15,200);
alpha	   	= 1e-15;
cgit       	= 8; 
pcgTol     	= 1e-1;
maxIter    	= 10;
minUpdate 	= 1e-2;
Iact 		= sparse(I,prod(n),prod(n));
boundsLow  	= (minimum(mtrue)-0.05)*ones(size(Iact,2)); # minimum(mtrue) is sea...
boundsHigh 	= (maximum(mtrue)+0.05)*ones(size(Iact,2));
maxStep		= 0.1*maximum(boundsHigh);

mtrue   	= Iact'*mtrue[:];
mref    	= Iact'*mref[:];
mback       = mref - Iact*(Iact'*mref)

regparams = [1.0,1.0,1.0,1e-6];
regfun(m,mref,Minv) = wdiffusionReg(m,mref,Minv,Iact=Iact,C = regparams);





if invertJoint
	filenamePrefix = "SEG_Joint_"
elseif invertEik
	filenamePrefix = "SEG_Eik_"
elseif invertDC
	filenamePrefix = "SEG_DC_"
end


function dumpResults(mc,Dc,iter,pInv,pMis)
	mDc = reshape(Iact*vec(mc)+mback,tuple(Minv.n...))	
	if plotting
		figure(888,figsize = (22,10));
		clf();
	end
	file = string(filenamePrefix,Minv.n,iter,".dat");
	if plotting 
		plotModel(mDc,false,[],0,[minimum(boundsLow),maximum(boundsHigh)],file);
	end
	writedlm(file,convert(Array{Float16},mDc[:]));
end
dumpResults(mref,[],"ref",[],[]);
dumpResults(mtrue,[],"true",[],[]);




###########################################################################################################################
###########################################################################################################################
################ DC Resistivity Alone Code ################################################################################
###########################################################################################################################


if invertDC
	Mfwd    = getRegularMesh(MinvNodal.domain,nfwd);
	sigtrue = velToConductMod(mtrue,3.0,0.1,1.0)[1]
	sigref  = velToConductMod(mref,3.0,0.1,1.0)[1]
	sigback      = sigref - Iact*(Iact'*sigref)
	
	@everywhere begin
		modfunDC(m) = velToConductMod(m,3.0,0.1,1.0);
	end
	

	println("------- setup pFor -----")
	Ainv = getMUMPSsolver()
	Mesh2Mesh  = getInterpolationMatrix(Minv,Mfwd)'
	redoSetup = true; # load setup if possible
	redoData =  true; # if possible, do not recompute data

	println("---------  sending out pMis ---------- ")
	
	
	if nworkers()>1
		# Julia sets BLAS #threads to be 1 if workers are used.
		set_num_threads(nworkers());
	end
	pFor,Sources,Receivers,Wd,dobsDC,dobs0 = setupDivSigGrad(redoSetup,redoData,Mfwd,Ainv,Mesh2Mesh,sigtrue,sigref);
	Wt   = 1./(mean(abs.(vec(dobsDC)))/2+abs.(dobsDC));
	dobsDC += 0.01*randn(size(dobsDC))*mean(abs.(vec(dobsDC)))
	println("adding noise of magnitude: ",0.01*randn()*mean(abs.(vec(dobsDC))), " Compared to ",mean(abs.(vec(dobsDC))));
	gloc = GlobalToLocal(Iact'*Mesh2Mesh,Mesh2Mesh'*sigback);
		
	
	if invertJoint
		worker = workers()[1];
		pMis = initRemoteChannel(getMisfitParam,worker, pFor,Wt,dobsDC,SSDFun,modfunDC,gloc);
		pMis = [pMis];
		modfun = identityMod;
	elseif invertDC
		pMis = getMisfitParam(pFor,Wt,dobsDC,SSDFun,identityMod,gloc);
		modfun = modfunDC;
	end
	pInv = getInverseParam(Minv,modfun,
                regfun,alpha,mref[:],
                boundsLow,boundsHigh,
                maxStep=maxStep,pcgMaxIter=cgit,pcgTol=pcgTol,
                minUpdate=minUpdate, maxIter = maxIter,HesPrec=HesPrec)
	println("-------- solve DivSigGrad ---------")	
end
###########################################################################################################################
###########################################################################################################################
################ Travel Time Tomography Code ##############################################################################
###########################################################################################################################

if invertEik
	jump 	= 10; 
	offset 	= 128;
	pad 	= 2;
	
	slowref 	= velocityToSlowSquared(mref)[1];
	slowback	= slowref - Iact*(Iact'*slowref);
	
	dataFilenamePrefix = string("DATA_SEG",tuple((MinvNodal.n+1)...));
	resultsFilename    = string("travelTimeInvSEG",tuple((MinvNodal.n+1)...),".dat");
	prepareTravelTimeDataFiles(mtrue,MinvNodal,mref,boundsHigh,boundsLow,dataFilenamePrefix,pad,jump,offset);
	
	if invertJoint
		pMisEik = setupTravelTimeTomography(dataFilenamePrefix, MinvNodal,Iact,slowback,velocityToSlowSquared)[3];
		modfun = identityMod;
		pMis = [pMis;pMisEik];
	elseif invertEik
		pMis = setupTravelTimeTomography(dataFilenamePrefix, MinvNodal,Iact,slowback,identityMod)[3];
		modfun = velocityToSlowSquared;
	end	
	pInv = getInverseParam(Minv,modfun,regfun,alpha,mref[:],boundsLow,boundsHigh,
                     maxStep=maxStep,pcgMaxIter=cgit,pcgTol=pcgTol,
					 minUpdate=minUpdate, maxIter = maxIter,HesPrec=HesPrec);
end
#### Projected Gauss Newton
mc,Dc,flag = projGNCG(copy(mref[:]),pInv,pMis,dumpResults = dumpResults);
#### Barrier Gauss Newton
# mc,Dc,flag = barrierGNCG(copy(mref[:]),pInv,pMis,dumpResults = dumpResults,epsilon = 0.3);
mv("jInv.out",string(filenamePrefix,".out"),remove_destination=true);
t = 1.0;

