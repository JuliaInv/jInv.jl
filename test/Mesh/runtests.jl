@testset "Mesh" begin
	@testset "display" begin include("display.jl") end
	@testset "Constraints" begin include("testConstraints.jl") end
	@testset "Regular vs. Tensor" begin include("regularVStensor.jl") end
	@testset "testDiffOps.jl" begin include("testDiffOps.jl") end
	@testset "testInterpolationMatrix.jl" begin include("testInterpolationMatrix.jl") end
	@testset "testAvgOps.jl" begin include("testAvgOps.jl") end
	@testset "testBoundaryNodes.jl" begin include("testBoundaryNodes.jl") end
	@testset "testPersistency.jl" begin include("testPersistency.jl") end
	@testset "testIntPolyChain.jl" begin include("testIntPolyChain.jl") end
	@testset "testMassMatrices.jl" begin include("testMassMatrices.jl") end
	@testset "testForward.jl" begin include("testForward.jl") end
end
