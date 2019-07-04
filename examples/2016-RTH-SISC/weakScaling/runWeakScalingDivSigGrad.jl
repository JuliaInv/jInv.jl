using MAT
using jInv.Mesh
using DivSigGrad
using jInv.LinearSolvers
using jInv.ForwardShare
using jInv.InverseSolve
using jInv.Utils
using ArgParse
using KrylovMethods

include("getDCResistivitySourcesAndRecAll.jl")

function parse_commandline()
	s = ArgParseSettings()

	@add_arg_table s begin

	"--n1"
		help = "number of cells in x1 direction "
		arg_type = Int
		default = 48
	"--n2"
		help = "number of cells in x2 direction "
		arg_type = Int
		default = 48
	"--n3"
		help = "number of cells in x3 direction "
		arg_type = Int
		default = 24
	"--srcSpacing"
		help = "spacing between sources"
		arg_type = Int
		default = 2
	"--nthreads"
		help = "number of threads"
		arg_type = Int
		default = 1
	"--out"
		help ="file to write into"
		default = "times.csv"
	"--solver"
		help = "PDE solver: 1->PCG, 2->BlockPCG, 3->MUMPS"
		arg_type = Int
		default = 1
	end
	return parse_args(s)
end

function main()
	parsed_args = parse_commandline()

	nInv      = [128;128;64]
	m         = 1500 + 3000*rand(tuple(nInv...))
	# random conductivity model can be replaced by SEG model described in
	#
	# Aminzadeh, F., Brac, J., and Kunz, T., 1997. 3D Salt and Overthrust models. SEG/EAGE Modeling Series, No. 1: Distribution CD of Salt and
	# Overthrust models, SEG Book Series Tulsa, Oklahoma
	#
	# matfile   = matread("3Dseg12812864.mat")
	# m         = matfile["VELc"]


	n1 = parsed_args["n1"]
	n2 = parsed_args["n2"]
	n3 = parsed_args["n3"]
	ns = parsed_args["srcSpacing"]
	nthreads = parsed_args["nthreads"]

	out = parsed_args["out"]

	# set number of threads for openblac
	set_num_threads(nthreads)

	solver = parsed_args["solver"]

	domain  = [0;13.5;0;13.5;0;4.2]
	n       = [n1;n2;n3]

	M    = getRegularMesh(domain,n)
	MInv = getRegularMesh(domain,nInv)
	Q,P  = getDCResistivitySourcesAndRecAll(M,srcSpacing=[ns,ns])
	Q    = Q[:,1:10]
	nsrc = size(Q,2)

	@printf "problem size: n1=%d\tn2=%d\tn3=%d\tns=%d\n" n1 n2 n3 nsrc
	@printf "nworkers()=%d\n" nworkers()
	@printf "nthreads=%d \n" nthreads

	if solver==3
		@printf "solver=MUMPS\n\n"
		@printf "warmup\ttotal\tfactorization\tsolve\n"
	elseif solver==2
		@printf "solver=BlockPCG\n\n"
		@printf "warmup\ttotal\titer\tprecond\tcgTime\tmatvec\n"
	elseif solver==1
		@printf "solver=PCG\n\n"
		@printf "warmup\ttotal\titer\tprecond\tcgTime\tmatvec\n"
	end
	Ainv = []
	if solver==1
		Ainv   = getIterativeSolver(KrylovMethods.cg,out=-1,sym=1,PC=:jac);
		Ainv.tol=1e-8
		Ainv.nthreads=nthreads
	elseif solver==2
 		Ainv         = getBlockIterativeSolver(KrylovMethods.blockCG,out=-1,sym=1,PC=:jac);
		Ainv.maxIter = 1000
		Ainv.tol     = 1e-8
		Ainv.nthreads=nthreads
		Ainv.out     =-1
	elseif solver==3
		Ainv = getMUMPSsolver()
	elseif solver 		== 4
		levels      	= 5;
		numCores 		= nthreads;
		maxIter     	= 20;
		relativeTol 	= 1e-8;
		relaxType   	= "SPAI";
		relaxParam  	= 1.0;
		relaxPre 		= 2;
		relaxPost   	= 2;
		cycleType   	='V';
		coarseSolveType = "MUMPS";
		MG 				= getMGparam(levels,numCores,maxIter,relativeTol,relaxType,relaxParam,relaxPre,relaxPost,cycleType,coarseSolveType);
		Ainvt   			= getSA_AMGsolver(MG, "PCG",sym=1,out=-1);
		Ainv   			= getSA_AMGsolver(MG, "PCG",sym=1,out=-1);
	end

	pForp = Array{RemoteRef{Channel{Any}}}(nworkers())
	pFors = Array{RemoteRef{Channel{Any}}}(nworkers())
	Mesh2Mesh  = getInterpolationMatrix(MInv,M)'
	m          = Mesh2Mesh'*vec(m)
	glocp      = Array{RemoteRef{Channel{Any}}}(nworkers())

	for j=1:nworkers()
		pForp[j] = @spawnat workers()[j] DivSigGradParam(M,Q,P,[1.0],Ainv)
		pFors[j] = @spawnat workers()[j] DivSigGradParam(M,Q[:,1],P,[1.0],Ainv)
		glocp[j] = @spawnat workers()[j] sparse(I,M.nc,M.nc)
	end

	# warm up
	wTime = @elapsed begin
	res = getData(vec(m),pFors,glocp)
	end
	@printf "%3.2f\t" wTime
	clear!(res[1])
	clear!(res[2])
	clear!(pFors)

	tTime = @elapsed begin
	getData(vec(m),pForp,glocp)
	end
	@printf "%3.2f\t" tTime

	if solver==3
		fTime = 0.0; sTime =0.0
		for p=1:length(pForp)
			pf = take!(pForp[p])
			fTime += pf.Ainv.facTime
			sTime += pf.Ainv.solveTime
		end
		fTime /= nworkers()
		sTime /= nworkers()
		@printf "%3.4f\t\t%3.4f\n" fTime sTime

		f = open(out,"a")
		str = @sprintf "%d,%d,%d,%d,%d,%d,,%3.5f,%3.5f,%3.5f,%3.5f\n" nworkers() solver n1 n2 n3 nsrc wTime tTime fTime sTime
		write(f,str )
		close(f)
	elseif solver==1 || solver==2
		nIter = 0; pcTime = 0.0; cgTime=0.0; mvTime=0.0
		for p=1:length(pForp)
			pf = take!(pForp[p])
			nIter += pf.Ainv.nIter
			pcTime += pf.Ainv.timePC
			cgTime += pf.Ainv.timeSolve
			mvTime += pf.Ainv.timeMV
			clear!(pf)
		end
		nIter /= nworkers()
		pcTime /=nworkers()
		cgTime /=nworkers()
		mvTime /=nworkers()
		@printf "%5d\t%3.4f\t%3.4f\t%3.4f\n" nIter  pcTime cgTime mvTime

		f = open(out,"a")
		str = @sprintf "%d,%d,%d,%d,%d,%d,%3.5f,%3.5f,%3.5f,%3.5f,%3.5f,%d\n" nworkers() solver n1 n2 n3 nsrc wTime tTime pcTime cgTime mvTime nIter
		write(f,str )
		close(f)
	elseif solver == 4
		nIter = 0; tsetup = 0.0; tSol = 0.0
		for p=1:length(pForp)
			pf = take!(pForp[p])
			nIter  += pf.Ainv.nIter
			tsetup += pf.Ainv.timeSetup
			tSol   += pf.Ainv.timeSolve
			clear!(pf)
		end
		@printf "%5d\t%3.4f\t%3.4f\t\n" nIter tsetup tSol
	end

end
main()
