
@testset "ForwardShare" begin
    @testset "testGetDataParallel.jl" begin
        include("testGetDataParallel.jl")
    end
	@testset "testPrepareMesh2Mesh.jl" begin
        include("testPrepareMesh2Mesh.jl")
    end
    @testset "testTest.jl" begin
        include("testTest.jl")
    end
    @testset "testGetSensMat.jl" begin
        include("testGetSensMat.jl")
    end
    @testset "testClear.jl" begin
        include("testClear.jl")
    end
end
