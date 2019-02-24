@everywhere using FWI
using Test
# using PyPlot

include("setupFWItests.jl")
# include("../src/plotResults.jl")

Dobs, = getData(m,pFor)
# setup inversion parameters
Wd    = 1./(abs(Dobs) + 0.1)
alpha = 1e-3
mref  = ones(nx*nz)
pInv  = InversionParam(alpha,gamma,mref,Mesh)



opti = (GNinv, NLCGinv, SDinv)

println("Test optimizers\n")
for k=1:length(opti)
    println("\ttesting $(string(opti[k]))")
    
    mc,his = opti[k](copy(mref),Dobs,Wd,pFor,pInv,maxIter=4,doPlot=false)
    mp,hisp = opti[k](copy(mref),Dobs,Wd,pForp,pInv,maxIter=4)
    
    @test norm(mc-mp)/norm(mc) < 1e-10
    @test norm(vec(his-hisp))/norm(vec(his)) < 1e-10
    @test size(his,1)==5
    @test all(his[5,:] .!= 0) 
    @test all(diff(his[:,3]).<0) # descent method
    @test all(diff(his[:,4]).<0) # reduction in gradient
    
    println("\t\tOK!")
end
