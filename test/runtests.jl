using Base.Test

include("setupTests.jl")

@testset "jInv" begin
	include("Mesh/runtests.jl")
	include("ForwardShare/runtests.jl")
	include("Utils/runTests.jl")
	include("InverseSolve/runtests.jl")
	include("LinearSolvers/runtests.jl")
end
