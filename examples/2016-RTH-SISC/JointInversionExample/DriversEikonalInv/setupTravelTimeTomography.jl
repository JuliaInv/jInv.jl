function setupTravelTimeTomography(filenamePrefix::String, Minv,Iact,sback,modFun)
###########################################################
### Read receivers and sources files
RCVfile = string(filenamePrefix,"_rcvMap.dat");
SRCfile = string(filenamePrefix,"_srcMap.dat");

srcNodeMap = readSrcRcvLocationFile(SRCfile,Minv);
rcvNodeMap = readSrcRcvLocationFile(RCVfile,Minv);

Q = generateSrcRcvProjOperators(Minv.n+1,srcNodeMap);
Q = Q.*1/(norm(Minv.h)^2);
P = generateSrcRcvProjOperators(Minv.n+1,rcvNodeMap);

println("Travel time tomography: ",size(Q,2)," sources.");
#############################################################################################

HO = false;

println("Reading data:");

(DobsEik,WdEik) =  readDataFileToDataMat(string(filenamePrefix,"_travelTime.dat"),srcNodeMap,rcvNodeMap);

N = prod(Minv.n+1);

########################################################################################################
##### Set up remote workers ############################################################################
########################################################################################################

EikMPIWorkers = nworkers(); # this just set the maximal MPI workers. To activate parallelism, run addprocs()

(pFor,contDiv,SourcesSubInd) = getEikonalInvParam(Minv,Q,P,HO,EikMPIWorkers);


misfun = SSDFun

Wd 		= Array{Array{Float64}}(length(pFor));
dobs 	= Array{Array{Float64}}(length(pFor));
for i=1:length(pFor)
	I_i = SourcesSubInd[i];
	Wd[i]   = 0.5*WdEik[:,I_i];
	dobs[i] = DobsEik[:,I_i];
end
WdEik = 0;
DobsEik = 0;

pMisRFs = getMisfitParam(pFor, Wd, dobs, misfun, Iact,sback,ones(length(pFor)),modFun);

return (Q,P,pMisRFs,SourcesSubInd,contDiv);
end


