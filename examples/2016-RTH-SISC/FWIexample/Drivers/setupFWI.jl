function setupFWI(m,filenamePrefix::String,resultsOutputFolderAndPrefix::String,plotting::Bool,
		workersFWI::Array{Int64,1}=workers(),maxBatchSize::Int64 = 48,
		Ainv::AbstractSolver = getMUMPSsolver([],0,0,2), misfun::Function=SSDFun,useFilesForFields::Bool = false)

		
file = matopen(string(filenamePrefix,"_PARAM.mat"));
n_cells = read(file,"n");
OmegaDomain = read(file,"domain");
Minv = getRegularMesh(OmegaDomain,n_cells);
gamma = read(file,"gamma");
omega = read(file,"omega");
waveCoef = read(file,"waveCoef");
if length(omega)==1 # matlab saves a 1 variable array as scalar.
	omega = [omega];
	waveCoef = [waveCoef];
end

boundsLow = read(file,"boundsLow");
boundsHigh = read(file,"boundsHigh");
mref =  read(file,"mref");
close(file);

resultsFilename = string(resultsOutputFolderAndPrefix,tuple((Minv.n+1)...),".dat");


### Read receivers and sources files
RCVfile = string(filenamePrefix,"_rcvMap.dat");
SRCfile = string(filenamePrefix,"_srcMap.dat");

srcNodeMap = readSrcRcvLocationFile(SRCfile,Minv);
rcvNodeMap = readSrcRcvLocationFile(RCVfile,Minv);

Q = generateSrcRcvProjOperators(Minv.n+1,srcNodeMap);
Q = Q.*1/(norm(Minv.h)^2);
println("We have ",size(Q,2)," sources");
P = generateSrcRcvProjOperators(Minv.n+1,rcvNodeMap);


########################################################################################################
##### Set up remote workers ############################################################################
########################################################################################################

N = prod(Minv.n+1);

Iact = sparse(1.0I,N,N);
mback   = zeros(Float64,N);
## Setting the sea constant:
mask = zeros(N);
mask[abs.(m[:] .- minimum(m)) .< 1e-2] = 1.0;
mask[gamma[:] .>= 0.95*maximum(gamma)] = 1.0;
# setup active cells
mback = vec(mref[:].*mask);
sback = velocityToSlowSquared(mback)[1];
sback[mask .== 0.0] = 0.0;
Iact = Iact[:,mask .== 0.0];

boundsLow = Iact'*boundsLow;
boundsHigh = Iact'*boundsHigh;
mref = Iact'*mref[:];

misfun = SSDFun;

####################################################################################################################
####################################################################################################################

println("Reading FWI data:");

batch = min(size(Q,2),maxBatchSize);
(pForFWI,contDivFWI,SourcesSubIndFWI) = getFWIparam(omega,waveCoef,vec(gamma),Q,P,Minv,Ainv,workersFWI,batch,useFilesForFields);

# write data to remote workers
Wd   = Array{Array{Complex128,2}}(length(pForFWI))
dobs = Array{Array{Complex128,2}}(length(pForFWI))
for k = 1:length(omega)
	omRound = string(round((omega[k]/(2*pi))*100.0)/100.0);
	(DobsFWIwk,WdFWIwk) =  readDataFileToDataMat(string(filenamePrefix,"_freq",omRound,".dat"),srcNodeMap,rcvNodeMap);
	for i = contDivFWI[k]:contDivFWI[k+1]-1
		I_i = SourcesSubIndFWI[i]; # subset of sources for ith worker.
		Wd[i] 	= WdFWIwk[:,I_i];
		dobs[i] = DobsFWIwk[:,I_i];
	end
	DobsFWIwk = 0;
	WdFWIwk = 0;
end


pMisFWIRFs = getMisfitParam(pForFWI, Wd, dobs, misfun, Iact,sback);

########################################################################################################
##### Set up remote workers ############################################################################
########################################################################################################

SourcesSubInd    = SourcesSubIndFWI;
pMis 			 = pMisFWIRFs;
contDiv			 = contDivFWI;

return Q,P,pMis,SourcesSubInd,contDiv,Iact,sback,mref,boundsHigh,boundsLow,resultsFilename;
end
