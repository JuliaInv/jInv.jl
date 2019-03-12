using Test
using LinearAlgebra

include("setupTests.jl")

@testset "jInv" begin
	include("Mesh/runtests.jl")
	include("ForwardShare/runtests.jl")
	include("Utils/runtests.jl")
	include("InverseSolve/runtests.jl")
	include("LinearSolvers/runtests.jl")
end
