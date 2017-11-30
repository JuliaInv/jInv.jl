function prepareTravelTimeDataFiles(m,Minv::RegularMesh,mref,boundsHigh,boundsLow,filenamePrefix::String,pad::Int64,jump::Int64,offset::Int64,useFilesForFields::Bool = false)
### Here we generate data files for sources and receivers. This can be replaced with real sources/receivers files.
########################## m is in Velocity here. ###################################

RCVfile = string(filenamePrefix,"_rcvMap.dat");
SRCfile = string(filenamePrefix,"_srcMap.dat");
writeSrcRcvLocFile(SRCfile,Minv,pad,jump);
writeSrcRcvLocFile(RCVfile,Minv,pad,1);

dataFullFilename = string(filenamePrefix,"_travelTime");
HO = false;
prepareTravelTimeDataFiles(m, Minv, filenamePrefix,dataFullFilename,offset,HO,useFilesForFields);


file = matopen(string(filenamePrefix,"_PARAM.mat"), "w");
write(file,"MinvOmega",Minv.domain);
write(file,"boundsLow",boundsLow);
write(file,"boundsHigh",boundsHigh);
write(file,"mref",mref);
write(file,"MinvN",Minv.n);
write(file,"HO",HO);
close(file);
return;
end

function prepareTravelTimeDataFiles(m, Minv::RegularMesh, filenamePrefix::String,dataFullFilename::String, offset::Int64,HO::Bool,useFilesForFields::Bool = false)

########################## m is in Velocity here. ###################################

RCVfile = string(filenamePrefix,"_rcvMap.dat");
SRCfile = string(filenamePrefix,"_srcMap.dat");
srcNodeMap = readSrcRcvLocationFile(SRCfile,Minv);
rcvNodeMap = readSrcRcvLocationFile(RCVfile,Minv);

Q = generateSrcRcvProjOperators(Minv.n+1,srcNodeMap);
Q = Q.*(1/(norm(Minv.h)^2));
P = generateSrcRcvProjOperators(Minv.n+1,rcvNodeMap);

# compute observed data

println("~~~~~~~ Getting data Eikonal: ~~~~~~~");
(pForEIK,contDivEIK,SourcesSubIndEIK) = getEikonalInvParam(Minv,Q,P,HO,nworkers(),useFilesForFields);

(D,pForEIK) = getData(velocityToSlowSquared(m[:])[1],pForEIK,ones(length(pForEIK)),true);


Dobs = Array{Array{Float64,2}}(length(pForEIK))
for k = 1:length(pForEIK)
	Dobs[k] = fetch(D[k]);
end
Dobs = arrangeRemoteCallDataIntoLocalData(Dobs);

# D should be of length 1 becasue constMUSTBeOne = 1;
Dobs += 0.01*mean(abs.(Dobs))*randn(size(Dobs,1),size(Dobs,2));
Wd = (1.0./(abs.(Dobs)+ 0.1*mean(abs.(Dobs))));
Wd = limitDataToOffset(Wd,srcNodeMap,rcvNodeMap,offset);
writeDataFile(string(dataFullFilename,".dat"),Dobs,Wd,srcNodeMap,rcvNodeMap);

return (Dobs,Wd)
end