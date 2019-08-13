module LinearSolvers

	abstract type AbstractSolver end
	abstract type AbstractDirectSolver <: AbstractSolver end
	function factorLinearSystem! end
	export AbstractSolver, AbstractDirectSolver, factorLinearSystem!

	using KrylovMethods
	using Distributed
	using SparseArrays
	using LinearAlgebra
	using Pkg
	# check if ParSPMatVec is available
	global hasParSpMatVec = false
	try
		using ParSpMatVec;
		global hasParSpMatVec = ParSpMatVec.isBuilt();
	catch
	end


	# check if Pardiso is installed
	const minPardisoVersion = VersionNumber(0,1,2)
	global hasPardiso = false
	try
		using Pardiso
		global hasPardiso = true
		if myid()==1
		  vPardiso = Pkg.installed()["Pardiso"]
		  if vPardiso < minPardisoVersion
		    warn("jInv Pardiso support requires Pardiso.jl version $(minPardisoVersion) or greater. Pardiso support will not be loaded")
		    hasPardiso = false
		  end
		end
	catch
	end



	export solveLinearSystem!,solveLinearSystem

	solveLinearSystem(A,B,param::AbstractSolver,doTranspose::Int=0) = solveLinearSystem!(A,B,zeros(eltype(B),size(B)),param,doTranspose)

	import Distributed.clear!
	function clear!(M::AbstractSolver)
		M.Ainv = []
	end

	include("iterativeWrapper.jl")
	include("blockIterativeWrapper.jl")
	include("PardisoWrapper.jl")
	include("juliaWrapper.jl")

	export clear!

end # module LinearSolvers
