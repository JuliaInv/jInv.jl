module Mesh

using LinearAlgebra
using LinearAlgebra.BLAS
using jInv.Utils
using SparseArrays
using Printf
using Distributed

export AbstractMesh
export AbstractTensorMesh
abstract type AbstractMesh end
abstract type AbstractTensorMesh <: AbstractMesh end

import Distributed.clear!
function clear!(M::AbstractTensorMesh)
	M.Div  = clear!(M.Div )
	M.Grad = clear!(M.Grad)
	M.Curl = clear!(M.Curl)
	M.Af   = clear!(M.Af  )
	M.Ae   = clear!(M.Ae  )
	M.An   = clear!(M.An  )
	M.V    = clear!(M.V   )
	M.F    = clear!(M.F   )
	M.L    = clear!(M.L   )
	M.Vi   = clear!(M.Vi  )
	M.Fi   = clear!(M.Fi  )
	M.Li   = clear!(M.Li  )
	M.nLap = clear!(M.nLap  )
end

function speye(n)
	return sparse(1.0I,n,n);
end

function spdiagm(x::Vector)
	return sparse(Diagonal(x));
end


include("generic.jl")
include("tensor.jl")
include("regular.jl")
include("interpmat.jl")
include("display.jl")
include("regularVectorFaces.jl")

export getNodalConstraints, getEdgeConstraints, getFaceConstraints
getNodalConstraints(M::AbstractMesh) = (UniformScaling(1.0), UniformScaling(1.0), 1:prod(M.n+1))
getEdgeConstraints(M::AbstractMesh) = (UniformScaling(1.0), UniformScaling(1.0), 1:sum(M.ne))
getFaceConstraints(M::AbstractMesh) = (UniformScaling(1.0), UniformScaling(1.0))

export clear!

end
