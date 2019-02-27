
using Test
using jInv.ForwardShare
using jInv.LinearSolvers
using jInv.Mesh

mutable struct TestForwardProbType <: ForwardProbType
    M
    Sources
    Obs
    Fields
    Ainv
end


msh = getTensorMesh3D(ones(3), ones(3), ones(3))

Ainv = getJuliaSolver()

pf = TestForwardProbType(msh,ones(3),ones(3),ones(3),Ainv);

clear!(pf; clearAll=true)

@test isempty(pf.M)
@test isempty(pf.Sources)
@test isempty(pf.Obs)
@test isempty(pf.Fields)