export GlobalToLocal, getGlobalToLocal
export interpLocalToGlobal, interpGlobalToLocal
export prepareGlobalToLocal

"""
mutable struct jInv.InverseSolve.GlobalToLocal

Maps global model to local model

sigLocal = PForInv*sigGlobal + sigmaBackground

Fields:
	PForInv::SparseMatrixCSC - linear operator between inverse mesh and
	                           forward mesh, e.g., linear interpolation matrix
                               and projection on active set.
	sigmaBackground::Vector  - background model

Constructors:
    getGlobalToLocal(P::SparseMatrixCSC)
	getGlobalToLocal(P::SparseMatrixCSC,sigmaBack::Vector)

Example:
	Mesh2Mesh = getInterpolationMatrix(Minv,Mfwd)   # Minv and Mfwd have different resolutions
	sigmaBack = 1.2*ones(Minv.nc)                   # put background conductivity
	gloc      = getGlobalToLocal(Mesh2Mesh,sigmaBack)

	sigLocal     = gloc.PForInv\' * sigGlobal + gloc.sigmaBack
	# this call is equivalent, but needed in case the Mesh2Mesh matrix is compressed
	sigLocalFast = interpGlobalToLocal(sigGlobal,gloc.PForInv,gloc.sigmaBack)
"""
struct GlobalToLocal
	PForInv::Union{SparseMatrixCSC,AbstractFloat,AbstractModelTransform} # interpolation matrix from fwd mesh to inv mesh
	sigmaBackground::Union{Vector{Float64},AbstractFloat,AbstractModel} #  (# of cells fwd mesh)
end # mutable struct GlobalToLocal

# Constructors
getGlobalToLocal(P) = GlobalToLocal(P,1e-8)
getGlobalToLocal(P,sigBack) = GlobalToLocal(P,sigBack)
getGlobalToLocal(P,sigBack,fname) = GlobalToLocal(P,sigBack)

function prepareGlobalToLocal(Mesh2Mesh,Iact,sigmaBackground,fname)
	return getGlobalToLocal(Iact'*Mesh2Mesh,interpGlobalToLocal(sigmaBackground,Mesh2Mesh),fname)
end
