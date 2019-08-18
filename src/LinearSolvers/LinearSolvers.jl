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


	export solveLinearSystem!,solveLinearSystem

	solveLinearSystem(A,B,param::AbstractSolver,doTranspose::Int=0) = solveLinearSystem!(A,B,zeros(eltype(B),size(B)),param,doTranspose)

	import Distributed.clear!
	function clear!(M::AbstractSolver)
		M.Ainv = []
	end

	include("iterativeWrapper.jl")
	include("blockIterativeWrapper.jl")
	include("juliaWrapper.jl")

	export clear!

end # module LinearSolvers
