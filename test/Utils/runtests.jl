
@testset "Utils" begin
    @testset "testSparseUtils.jl" begin
        include("testSparseUtils.jl")
    end
    @testset "testSortpermFast.jl" begin
        include("testSortpermFast.jl")
    end
    @testset "testUniqueIdx.jl" begin
        include("testUniqueIdx.jl")
    end
    @testset "testVariousUtils.jl" begin
        include("testVariousUtils.jl")
    end
    @testset "testTests.jl" begin
        include("testTests.jl")
    end
end
