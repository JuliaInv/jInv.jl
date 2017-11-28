function prepareFWIDataFiles(m,Minv::RegularMesh,mref,boundsHigh,boundsLow,
							filenamePrefix::String,omega::Array{Float64,1},waveCoef::Array{Complex128,1},
							pad::Int64,ABLpad::Int64,jump::Int64,offset::Int64=prod(Minv.n+1),workerList = workers(), 
							maxBatchSize::Int64=48, Ainv::AbstractSolver = getMUMPSsolver([],0,0,2),useFilesForFields::Bool = false)

########################## m is in Velocity here. ###################################
RCVfile = string(filenamePrefix,"_rcvMap.dat");
SRCfile = string(filenamePrefix,"_srcMap.dat");
writeSrcRcvLocFile(SRCfile,Minv,ABLpad,jump);
writeSrcRcvLocFile(RCVfile,Minv,ABLpad,1);
	
dataFullFilenamePrefix = string(filenamePrefix,"_freq");
gamma = prepareFWIDataFiles2(m, Minv, filenamePrefix,dataFullFilenamePrefix,omega,waveCoef,pad,ABLpad,offset,workerList,maxBatchSize,Ainv);	

file = matopen(string(filenamePrefix,"_PARAM.mat"), "w");
write(file,"boundsLow",boundsLow);
write(file,"boundsHigh",boundsHigh);
write(file,"mref",mref);
write(file,"domain",Minv.domain);
write(file,"n",Minv.n);
write(file,"gamma",gamma);
write(file,"pad",pad);
write(file,"omega",omega);
write(file,"waveCoef",waveCoef);
close(file);
return gamma;	
end


function prepareJointTravelTimeAndFWIDataFiles(m,Minv::RegularMesh,mref,boundsHigh,boundsLow,
							filenamePrefix::String,omega::Array{Float64,1},waveCoef::Array{Complex128,1},
							pad::Int64,ABLpad::Int64,jump::Int64,offset::Int64=prod(Minv.n+1),workerList = workers(), 
							maxBatchSize::Int64=48, Ainv::AbstractSolver = getMUMPSsolver([],0,0,2),useFilesForFields::Bool = false)

########################## m is in Velocity here. ###################################
RCVfile = string(filenamePrefix,"_rcvMap.dat");
SRCfile = string(filenamePrefix,"_srcMap.dat");
writeSrcRcvLocFile(SRCfile,Minv,ABLpad,jump);
writeSrcRcvLocFile(RCVfile,Minv,ABLpad,1);
	
dataFullFilenamePrefix = string(filenamePrefix,"_freq");
gamma = prepareFWIDataFiles2(m, Minv, filenamePrefix,dataFullFilenamePrefix,omega,waveCoef,pad,ABLpad,offset,workerList,maxBatchSize,Ainv,useFilesForFields);	
dataFullFilenamePrefix = string(filenamePrefix,"_travelTime");
HO = false;
prepareTravelTimeDataFiles(m, Minv, filenamePrefix,dataFullFilenamePrefix,offset,HO,useFilesForFields);

file = matopen(string(filenamePrefix,"_PARAM.mat"), "w");
write(file,"boundsLow",boundsLow);
write(file,"boundsHigh",boundsHigh);
write(file,"mref",mref);
write(file,"MinvOmega",Minv.domain);
write(file,"MinvN",Minv.n);
write(file,"gamma",gamma);
write(file,"pad",pad);
write(file,"omega",omega);
write(file,"waveCoef",waveCoef);
write(file,"HO",HO);
close(file);			
end


function prepareFWIDataFiles2(m, Minv::RegularMesh, filenamePrefix::String,dataFullFilenamePrefix::String, omega::Array{Float64,1}, 
								waveCoef::Array{Complex128,1}, pad::Int64,ABLpad::Int64,offset::Int64,workerList::Array{Int64,1},maxBatchSize::Int64,
								Ainv::AbstractSolver,useFilesForFields::Bool = false)


println("maxOmega*maxKappaSq*h: should be below 0.6");
println(omega[end]*maximum(Minv.h)*(1./minimum(m)));


RCVfile = string(filenamePrefix,"_rcvMap.dat");
SRCfile = string(filenamePrefix,"_srcMap.dat");
srcNodeMap = readSrcRcvLocationFile(SRCfile,Minv);
rcvNodeMap = readSrcRcvLocationFile(RCVfile,Minv);

Q = generateSrcRcvProjOperators(Minv.n+1,srcNodeMap);
P = generateSrcRcvProjOperators(Minv.n+1,rcvNodeMap);
Q = Q.*1/(norm(Minv.h)^2);
println("We have ",size(Q,2)," sources");
# compute observed data

ABLamp = getMaximalFrequency(1./(minimum(m).^2),Minv);
gamma = getABL(Minv,true,ones(Int64,Minv.dim)*ABLpad,ABLamp);
attenuation = 0.01*ABLamp;
gamma += attenuation; # adding Attenuation.

println("~~~~~~~ Getting data FWI: ~~~~~~~");
	
# Solve forward problem (should be relaced in reading data from file)

batch = min(size(Q,2),maxBatchSize);
(pFor,contDiv,SourcesSubInd) = getFWIparam(omega,waveCoef,vec(gamma),Q,P,Minv,Ainv,workerList,batch,useFilesForFields);

(D,pFor) = getData(velocityToSlowSquared(m[:])[1],pFor,ones(length(pFor)),true);

nsrc = size(Q,2);
nrcv = size(P,2);


for k = 1:length(omega)
	I = contDiv[k] : contDiv[k+1] - 1 
	Dobsk = Array{Array{Complex128,2}}(length(I));
	for i = 1:length(I)
		Dobsk[i] = fetch(D[I[i]]);
	end
	Dobsk = arrangeRemoteCallDataIntoLocalData(Dobsk);
	omRound = string(round((omega[k]/(2*pi))*100.0)/100.0);
	Wd_k = (1./(abs.(real(Dobsk))+0.1*mean(abs.(Dobsk)))) + 1im*(1./(abs.(imag(Dobsk))+0.1*mean(abs.(Dobsk))));
	# Wd_k = (1./(0.0*abs(real(Dobsk))+1.0*mean(abs(Dobsk)))) + 1im*(1./(0.0*abs(imag(Dobsk))+1.0*mean(abs(Dobsk))));
	Wd_k = limitDataToOffset(Wd_k,srcNodeMap,rcvNodeMap,offset);
	filename = string(dataFullFilenamePrefix,omRound,".dat");
	writeDataFile(filename,Dobsk,Wd_k,srcNodeMap,rcvNodeMap);
end
return gamma;
end