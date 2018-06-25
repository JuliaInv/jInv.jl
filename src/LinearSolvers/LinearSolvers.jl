module LinearSolvers

	abstract type AbstractSolver end
	abstract type AbstractDirectSolver <: AbstractSolver end
	function factorLinearSystem! end
	export AbstractSolver, AbstractDirectSolver, factorLinearSystem!

	using KrylovMethods

	# check if MUMPS can be used
	const minMUMPSversion = VersionNumber(0,0,1)
	hasMUMPS=false
	try
		using MUMPS
		hasMUMPS = true
		if myid()==1
			vMUMPS = Pkg.installed("MUMPS")
			if vMUMPS < minMUMPSversion;
				warn("MUMPS is outdated! Please update to version $(minMUMPSversion) or higher.")
			end
		end
	catch
	end

	# check if Pardiso is installed
	const minPardisoVersion = VersionNumber(0,1,2)
	hasPardiso = false
	try
		using Pardiso
		hasPardiso = true
		if myid()==1
		  vPardiso = Pkg.installed("Pardiso")
		  if vPardiso < minPardisoVersion
		    warn("jInv Pardiso support requires Pardiso.jl version $(minPardisoVersion) or greater. Pardiso support will not be loaded")
		    hasPardiso = false
		  end
		end
	catch
	end

	# check if ParSPMatVec is available
	hasParSpMatVec = false
	const minVerParSpMatVec = VersionNumber(0,0,1)
		vParSpMatVec = VersionNumber(0,0,0)
	try
 		using ParSpMatVec
 		hasParSpMatVec = true
 		if myid()==1
 			vParSpMatVec = Pkg.installed("ParSpMatVec")
 			if vParSpMatVec < minVerParSpMatVec;
 				warn("ParSpMatVec is outdated! Please update to version $(minVerParSpMatVec) or higher.")
 			end
 		end
	catch
	end
	if hasParSpMatVec
		using ParSpMatVec
	end

	export solveLinearSystem!,solveLinearSystem

	solveLinearSystem(A,B,param::AbstractSolver,doTranspose::Int=0) = solveLinearSystem!(A,B,zeros(eltype(B),size(B)),param,doTranspose)

	import Base.clear!
	function clear!(M::AbstractSolver)
		M.Ainv = []
	end

	include("mumpsWrapper.jl")
	include("iterativeWrapper.jl")
	include("blockIterativeWrapper.jl")
	include("PardisoWrapper.jl")
	include("juliaWrapper.jl")

	export clear!

end # module LinearSolvers
