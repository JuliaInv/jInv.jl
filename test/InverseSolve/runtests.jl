
@testset "InverseSolve" begin
    @testset "testModels.jl" begin
        include("testModels.jl")
    end
    @testset "testMisfit.jl" begin
        include("testMisfit.jl")
    end
    @testset "testRegularizers.jl" begin
       include("testRegularizers.jl")
    end
    @testset "testMisfitParams.jl" begin
       include("testMisfitParams.jl")
    end
    @testset "testLeastSquares.jl" begin
       include("testLeastSquares.jl")
    end
    @testset "testGetHessian.jl" begin
       include("testGetHessian.jl")
    end
    @testset "testHessMatVec.jl" begin
        include("testHessMatVec.jl")
    end
end
