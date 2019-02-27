using  jInv.Mesh
using  jInv.Utils
using  jInv.LinearSolvers
using  jInv.InverseSolve
using  jInvVis
using  EikonalInv
using  MAT
using  FWI
using  ForwardHelmholtz
using  PyPlot
#############################################################################################################


####################################################################################################################################
####################################################################################################################################
#######################   Problem definition and data creation #####################################################################
####################################################################################################################################
####################################################################################################################################

include("ex2DFWI/readModelAndGenerateMeshMref.jl");
include("ex2DFWI/prepareFWIDataFiles.jl");

m = readdlm("model/SEGmodel2Dsalt.dat");
m = m*1e-3;
m = m';

newSize          = [256,128];
pad     	     = 10;
ABLPad 		     = pad + 8;
jumpSrc 	 	 = 5
maxBatchSize     = 256;
omega   	     = [0.5,0.75,1.25,1.75]*2*pi;

offset  = ceil(Int64,(newSize[1]*(10.0/13.5)));
println("Offset is: ",offset)
domain = [0.0,13.5,0.0,4.2];

(m,Minv,mref,boundsHigh,boundsLow) = readModelAndGenerateMeshMref(m,pad,newSize,domain);


	limits = [1.5,4.5];
	figure(1,figsize = (22,10))
	plotModel(m,true,Minv,pad,limits);
	figure(2,figsize = (22,10))# ,figsize = (22,10)
	plotModel(mref,true,Minv,pad,limits);


useFilesForFields = false;

# ###################################################################################################################
dataFilenamePrefix = "ex2DFWI/DATA_SEG2D";

#######################################################################################################################


######################## DIRECT SOLVER #################################################
numCores 	= 4;
BLAS.set_num_threads(numCores);
Ainv = getMUMPSsolver([],0,0,2); # Alternatives: Ainv = getJuliaSolver(); 

##########################################################################################

println("omega*maximum(h): ",omega*maximum(Minv.h)*sqrt(maximum(1./(boundsLow.^2))));

# This is a list of workers for FWI. Ideally they should be on different machines.
workersFWI = [workers()[1]];
prepareFWIDataFiles(m,Minv,mref,boundsHigh,boundsLow,dataFilenamePrefix,omega,one(Complex128)*ones(size(omega)), pad,ABLPad,jumpSrc,
					offset,workersFWI,maxBatchSize,Ainv,useFilesForFields);
		
####################################################################################################################################
####################################################################################################################################
####################### Inversion setup for a single machine  ######################################################################
####################################################################################################################################
####################################################################################################################################		
dataFilenamePrefix = "ex2DFWI/DATA_SEG2D";
file = matread(string(dataFilenamePrefix,"_PARAM.mat"));
n_cells 		= file["n"];
domain 	        = file["domain"];
gamma 			= file["gamma"];
omega 			= file["omega"];
waveCoef 		= file["waveCoef"];
boundsLow 		= file["boundsLow"];
boundsHigh 		= file["boundsHigh"];
mref 			= file["mref"];
file = 0;

Minv = getRegularMesh(domain,n_cells);
if length(omega)==1 # matlab saves a 1 variable array as scalar.
	omega = [omega];
	waveCoef = [waveCoef];
end

### Read receivers and sources files
RCVfile = string(dataFilenamePrefix,"_rcvMap.dat");
SRCfile = string(dataFilenamePrefix,"_srcMap.dat");
srcNodeMap = readSrcRcvLocationFile(SRCfile,Minv);
rcvNodeMap = readSrcRcvLocationFile(RCVfile,Minv);
Q = generateSrcRcvProjOperators(Minv.n+1,srcNodeMap); 
Q = Q.*1/(norm(Minv.h)^2);
println("We have ",size(Q,2)," sources");
P = generateSrcRcvProjOperators(Minv.n+1,rcvNodeMap);


### Read the data files to an array of array pointers
Wd   = Array{Array{Complex128,2}}(length(omega))
dobs = Array{Array{Complex128,2}}(length(omega))
for k = 1:length(omega)
	omRound = string(round((omega[k]/(2*pi))*100.0)/100.0);
	(DobsFWIwk,WdFWIwk) =  readDataFileToDataMat(string(dataFilenamePrefix,"_freq",omRound,".dat"),srcNodeMap,rcvNodeMap);
	Wd[k] 	= WdFWIwk;
	dobs[k] = DobsFWIwk;
end



# setup active cells and background for the inversion.
N = prod(Minv.n+1);
Iact = sparse(I,N,N);
mback   = zeros(Float64,N);

########################################################################################################
##### Set up remote workers ############################################################################
########################################################################################################

batch = size(Q,2);         # We solve all the sources at once.
useFilesForFields = false; # a flag whether to store the fields on the disk.			
## Set up a MUMPS direct solver #Ainv = getMUMPSsolver([],0,0,2); #
Ainv = getJuliaSolver(); 
## Choose the workers for FWI (here, its a single worker)
workersFWI = [1];
## Set up workers and division to tasks per frequencies ######################
(pFor,contDiv,SourcesSubInd) = getFWIparam(omega,waveCoef,vec(gamma),Q,P,Minv,Ainv,workersFWI,batch,useFilesForFields);
misfun = SSDFun;
pMis = getMisfitParam(pFor, Wd, dobs, misfun, Iact,mback);


########################################################################################################
# Set up Inversion #################################################################################
########################################################################################################



## models are usually given in velocity (km/sec). We invert here for the slowness, which is 1/v.
mref 				= velocityToSlow(mref)[1];
t    				= copy(boundsLow);
boundsLow 			= velocityToSlow(boundsHigh)[1];
boundsHigh 			= velocityToSlow(t)[1]; t = 0;
modfun 				= slowToSlowSquared;


maxStep				=0.05*maximum(boundsHigh);

regparams 			= [1.0,1.0,1.0,1e-5];
cgit 				= 8; 
alpha 				= 1e-10;
pcgTol 				= 1e-1;
maxit 				= 8;
HesPrec 			= getExactSolveRegularizationPreconditioner();
regfun(m,mref,M) 	= wdiffusionRegNodal(m,mref,M,Iact=Iact,C=regparams);


pInv = getInverseParam(Minv,modfun,regfun,alpha,mref[:],boundsLow,boundsHigh,
                         maxStep=maxStep,pcgMaxIter=cgit,pcgTol=pcgTol,
						 minUpdate=1e-3, maxIter = maxit,HesPrec=HesPrec);

						 
plotting = true;			 
function plotIntermediateResults(mc,Dc,iter,pInv,PMis,resultsFilename="")
	# Models are usually shown in velocity.
	fullMc = slowSquaredToVelocity(reshape(Iact*pInv.modelfun(mc)[1] + mback,tuple((pInv.MInv.n+1)...)))[1];
	if plotting
		close(888);
		figure(888);
		plotModel(fullMc,false,[],0,[1.5,4.8]);
	end
end
						 
# Run one sweep of a frequency continuation procedure.
mc,Dc,flag,His = freqCont(copy(mref[:]), pInv, pMis,contDiv, 4, "",plotIntermediateResults,"Joint",1,1,"projGN");

subplot(1,3,1)
semilogy(His[end].F,"-o")
xlabel("PGNCG iterations")
title("misfit")

subplot(1,3,2)
semilogy(His[end].Rc,"-o")
xlabel("PGNCG iterations")
title("regularizer");

subplot(1,3,3)
semilogy(His[end].dJ,"-o")
xlabel("PGNCG iterations")
title("norm of proj grad");


