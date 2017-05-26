@testset "Mesh" begin
	@testset "display" begin include("display.jl") end
	@testset "Constraints" begin include("testConstraints.jl") end
	@testset "Regular vs. Tensor" begin include("regularVStensor.jl") end
	@testset "tbd" begin include("testDiffOps.jl") end
	@testset "tbd" begin include("testInterpolationMatrix.jl") end
	@testset "tbd" begin include("testAvgOps.jl") end
	@testset "tbd" begin include("testBoundaryNodes.jl") end
	@testset "tbd fail" begin include("testPersistency.jl"); @test 3==5; end
	@testset "tbd" begin include("testIntPolyChain.jl") end 
	@testset "tbd" begin include("testMassMatrices.jl") end
end

