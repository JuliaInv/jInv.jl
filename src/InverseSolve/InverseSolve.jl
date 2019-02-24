module InverseSolve

using KrylovMethods
using jInv.Utils
using jInv.Mesh
using jInv.ForwardShare
using jInv.LinearSolvers

import jInv.ForwardShare.ForwardProbType

export getName, AbstractMisfit

abstract type  AbstractMisfit end

export AbstractModel,AbstractModelDerivative,AbstractModelTransform

abstract type AbstractModel end
abstract type AbstractModelDerivative end
abstract type AbstractModelTransform end

using Distributed
using SparseArrays
using Printf
using LinearAlgebra

include("HessianPreconditioners.jl")

include("misfitParam.jl")

"""
mutable struct jInv.InverseSolve.InverseParam

Type storing parameters for Inversion.

Fields:
    Minv::AbstractMesh
    modelfun              - model function (evaluated by main worker), see models.jl
    regularizer           - regularizer(s), see regularizer.jl
    alpha                 - regularization parameter(s)
    mref                  - reference model(s)
    boundsLow::Vector     - lower bounds for model
    boundsHigh::Vector    - upper bounds for model
    maxStep::Real         - maximum step in optimization
    pcgMaxIter::Int       - maximum number of PCG iterations
    pcgTol::Real          - tolerance for PCG
    minUpdate::Real       - stopping criteria
    maxIter::Int          - maximum number of iterations
    HesPrec               - preconditioner for the Hessian.
    metaData::Dict        - Optional dictonary for storeing additional information
Constructor:
    getInverseParam

Example:
    Minv = getRegularMesh(domain,n)
    modelfun = expMod
    regularizer(m,mref,Minv) = wdiffusionReg(m,mref,Minv)
    alpha   = 1e-3
    mref    = zeros(Minv.nc)
    pInv = getInverseParam(Minv,modelfun,regularizer,alpha,mref)
"""
mutable struct  InverseParam
    MInv::AbstractMesh
    modelfun::Function
    regularizer::Union{Function,Array{Function}}
    alpha::Union{Float64,Array{Float64}}
    mref::Array
    boundsLow::Vector
    boundsHigh::Vector
    maxStep::Real
    pcgMaxIter::Int
    pcgTol::Real
    minUpdate::Real
    maxIter::Int
    HesPrec::HessianPreconditioner
    metaData::Dict
end

function Base.display(pInv::InverseParam)
    println("---jInv.InverseSolve.InverseParam---")
    println("inverse mesh type:    $(typeof(pInv.MInv))")
    println("number of cells:      $(pInv.MInv.nc)")
    println("model function:       $(pInv.modelfun)")
    println("maxStep:              $(pInv.maxStep)")
    println("pcgMaxIter:           $(pInv.pcgMaxIter)")
    println("pcgTol:               $(pInv.pcgTol)")
    println("minUpdate:            $(pInv.minUpdate)")
    println("maxIter:              $(pInv.maxIter)")
end

"""
function jInv.InverseSolve.getInverseParam(...)

Constructs an InverseParam

Required Input:

    Minv::AbstractMesh    - mesh of model
    modFun::Function      - model
    regularizer::Function - regularizer, see regularizer.jl
    alpha::Real           - regularization parameter
    mref                  - reference model
    boundsLow::Vector     - lower bounds for model
    boundsHigh::Vector    - upper bounds for model

Optional Inputs:

    maxStep::Real=1.0     - maximum step in optimization
    pcgMaxIter::Int=10    - maximum number of PCG iterations
    pcgTol::Real          - tolerance for PCG
    minUpdate::Real=1e-4  - stopping criteria
    maxIter::Int=10       - maximum number of iterations
"""
function getInverseParam(MInv::AbstractMesh,
                         modFun::Function,
                         regularizer::Union{Function,Array{Function}},
                         alpha::Union{Float64,Array{Float64}},
                         mref::Array,
                         boundsLow::Vector,
                         boundsHigh::Vector;
                         maxStep::Real=1.0,
                         pcgMaxIter::Int=10,
                         pcgTol::Real=1e-1,
                         minUpdate::Real=1e-4,
                         maxIter::Int=10,
                         HesPrec::HessianPreconditioner=getSSORRegularizationPreconditioner(1.0,1e-15,10),
                         metaData::Dict=Dict())

    return InverseParam(MInv, modFun, regularizer, alpha, mref, boundsLow, boundsHigh, maxStep,
                        pcgMaxIter, pcgTol, minUpdate, maxIter, HesPrec, metaData)
end



export InverseParam, getInverseParam

include("globalToLocal.jl")
include("models.jl")
include("misfit.jl")
include("regularizers.jl")
include("GNsolve.jl")
include("GNhis.jl")
include("projGN.jl")
include("projSD.jl")
include("projPCG.jl")
include("barrierGNCG.jl")
include("computeMisfit.jl")
include("computeGradMisfit.jl")
include("HessMatVec.jl")
include("iteratedTikhonov.jl")
include("getHessian.jl")
end
