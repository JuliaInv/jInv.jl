using MAT
using jInv.Mesh
using EikonalInv
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

	"--nsources"
		help = "number of sources per worker"
		arg_type= Int
		default = 36
	"--out"
		help ="file to write into"
		default = "times.csv"
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
	out = parsed_args["out"]
	nsources = parsed_args["nsources"]

	# set number of threads for openbla to one for Eikonal
	@everywhere set_num_threads(1)

	domain  = [0;13.5;0;13.5;0;4.2]
	n       = [n1;n2;n3]
	M    = getRegularMesh(domain,n-1)

	# interpolate model to right resolution
	if any(M.n .!= nInv)
		MInv = getRegularMesh(domain,nInv)
		m   = getInterpolationMatrix(MInv,getRegularMesh(domain,n))*vec(m)
	end

	# get sources and receivers
	RCVfile = "DATA_SEG(64,64,32)_rcvMap.dat";
	SRCfile = "DATA_SEG(64,64,32)_srcMap.dat";
	srcNodeMap = readSrcRcvLocationFile(SRCfile,M);
	rcvNodeMap = readSrcRcvLocationFile(RCVfile,M);
	Q = generateSrcRcvProjOperators(M.n+1,srcNodeMap);
	Q = Q.*(1/(norm(M.h)^2));
	Q = Q[:,1:nsources]
	P = generateSrcRcvProjOperators(M.n+1,rcvNodeMap);
	nsrc = size(Q,2)

	@printf "problem size: n1=%d\tn2=%d\tn3=%d\tns=%d\n" n1 n2 n3 nsrc
	@printf "nworkers()=%d\n" nworkers()

	# setup pFor
	HO = false
	Qp = kron(Q,ones(1,nworkers())); # dulicate sources so that all workers solve same problems
	Qs = kron(Q[:,1],ones(1,nworkers()))
	(pForp,contDivEIK,SourcesSubIndEIK) = getEikonalInvParam(M,Qp,P,HO,nworkers())
	(pFors,contDivEIK,SourcesSubIndEIK) = getEikonalInvParam(M,Qs,P,HO,nworkers())


	# warm up
	wTime = @elapsed begin
	res = getData(vec(m),pFors,ones(length(pFors)),true)
	end
	@printf "%3.2f\t" wTime
	clear!(res[1])
	clear!(res[2])
	clear!(pFors)

	tTime = @elapsed begin
		m = vec(m);
	for k=1:5
		getData(m,pForp,ones(length(pForp)),true)
	end
	end
	tTime /= 5.0;
	@printf "%3.2f\n" tTime
	f = open(out,"a")
	str = @sprintf "%d,%d,%d,%d,%d,,%3.5f,%3.5f\n" nworkers() n1 n2 n3 nsrc wTime tTime
	write(f,str )
	close(f)
end
main()
