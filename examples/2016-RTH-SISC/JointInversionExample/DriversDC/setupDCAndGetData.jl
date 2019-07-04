function setupDivSigGrad(redoSetup,redoData,pFor,Ainv,Mesh2Mesh,sigtrue,sigref)
println("------- setup or load experiment from disk ---------")
srcSpacing = [2,2]
srcPad     = [1,1,1,1]
setupFile = "testDivSigGradSEGsetupMUMPS.mat"
if isfile(setupFile)
	mf = matread(setupFile)
	if any(srcSpacing .!= mf["srcSpacing"]) || any(srcPad .!= mf["srcPad"]) || any(Mfwd.n.!=mf["nfwd"])|| any(Minv.n.!=mf["n"])
		println("source spacing or padding has changed. need to setup experiment")
		redoSetup = true
	end
else
	redoSetup = true
end
if redoSetup
	Sources,Receivers,Wd = getDCResistivitySourcesAndRecAll(nfwd,Mfwd,srcSpacing=srcSpacing,srcPad=srcPad)
	results = Dict("sigtrue"=>sigtrue,"Sources"=>Sources,"Receivers"=>Receivers,"Wd"=>Wd,"srcSpacing"=>srcSpacing,
					"srcPad"=>srcPad,"n"=>Minv.n,"nfwd"=>Mfwd.n)
	matwrite(setupFile,results)
else
	mf        = matread(setupFile)
	Sources   = mf["Sources"]
	Receivers = mf["Receivers"]
	Wd        = mf["Wd"]
end

nsrc        = size(Sources,2)
println("number of sources $nsrc")


println("------- generate or load data ----------")
# if worker == 0
	pFor = DivSigGradParam(Mfwd,Sources,Receivers,[],spzeros(0,0),Ainv);
# else
	# pFor = remotecall(worker,DivSigGradParam,Mfwd,Sources,Receivers,[],Ainv);
	# pFor = [pFor]
# end
dataFile = "testDivSigGradSEGdataMUMPS.mat"
if isfile(dataFile)
	mf = matread(dataFile)
	if size(Sources)!=size(mf["Sources"]) || norm(Sources-mf["Sources"],Inf)>1e-10
		println("sources have changed. need to recompute data.")
		redoData = true;
	end
	if size(sigtrue)!=size(mf["sigtrue"]) || norm(sigtrue-mf["sigtrue"],Inf)>1e-8
		println("sigma changed. need to recompute data.")
		redoData = true;
	end
        if length(mf["dobs"])!=size(Receivers,2)*nsrc
		println("receivers have changed. need to recompute data.")
		redoData = true;
	end
end

if redoData || !isfile(dataFile)
	println("get DC data for true model")

	@time dobsDC,     = getData(Mesh2Mesh'*vec(sigtrue),pFor)
	println("get DC data for homogeneous model")

	@time dobs0,     = getData(Mesh2Mesh'*vec(sigref),pFor)
	# if worker!=0
		# dobs0 = fetch(dobs0[1]);
		# dobsDC = fetch(dobsDC[1]);
	# end
	results = Dict("sigtrue"=>sigtrue,"sigref"=>sigref,"dobs"=>dobsDC,"dobs0"=>dobs0,"Sources"=>Sources)
	matwrite(dataFile,results)
else
	println("Reading data from file")
	mf     = matread(dataFile)
	dobsDC = mf["dobs"]
	dobs0  = mf["dobs0"]
end



return pFor,Sources,Receivers,Wd,dobsDC,dobs0
end
