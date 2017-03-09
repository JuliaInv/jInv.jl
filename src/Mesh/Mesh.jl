module Mesh

using Compat
using Base.BLAS
using jInv.Utils

export AbstractMesh
export AbstractTensorMesh
abstract AbstractMesh
abstract AbstractTensorMesh <: AbstractMesh

import Base.clear!
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

include("generic.jl")
include("tensor.jl")
include("regular.jl")
include("interpmat.jl")
include("display.jl")

export getEdgeConstraints, getFaceConstraints
getEdgeConstraints(M::AbstractMesh) = (UniformScaling(1.0), UniformScaling(1.0), 1:sum(M.ne))
getFaceConstraints(M::AbstractMesh) = (UniformScaling(1.0), UniformScaling(1.0))

export clear!

end
