export interpLocalToGlobal, interpGlobalToLocal
# Interpolate a vector x from the local mesh to the global mesh

"""
function jInv.ForwardShare.interpLocalToGlobal(mLocal,P,const)
function jInv.ForwardShare.interpGlobalToLocal(mGlobal,P,const)

convert a local (or global) model to a global (or local) mesh so that
    mGlobal = P * mLocal + const
    mLocal = P' * mGlobal + const

Examples:
    P = prepareMesh2Mesh(Msmall, Mbig, true)
    modelLoc =  interpGlobalToLocal(modelInv,P) # uncompresses the indices and does P'*modelInv

    P = prepareMesh2Mesh(Msmall, Mbig, true)
    modelInv =  interpGlobalToLocal(modelLoc,P) # uncompresses the indices and does P'*modelInv

Inputs:
    mLocal  - a model defined on the local (PDE) mesh
    mGlobal - a model defined on the global (inverse) mesh
    const   - a scalar constant shift (if not provided: const = 0.0)

Outputs:
    mGlobal - a model defined on the global (inverse) mesh
    mLocal  - a model defined on the local (PDE) mesh
"""
function interpLocalToGlobal(x::Vector{Float64}, P::AbstractFloat,y0::AbstractFloat=0.0)
	return P * x .+ y0
end

function interpGlobalToLocal(x::Vector{Float64}, P::AbstractFloat,y0::AbstractFloat=0.0)
	return P * x .+ y0
end

function interpGlobalToLocal(x,P)
    return P'*x
end

function interpGlobalToLocal(x,P,y0)
    return P'*x + y0
end

function interpLocalToGlobal(x,P)
    return P*x
end

function interpLocalToGlobal(x::Vector{Float64}, P::SparseMatrixCSC)

	if (eltype(P.nzval) == Int16) || (eltype(P.nzval) == Int8)
   		nzv = P.nzval
    	rv  = P.rowval
      	y   = zeros(P.m)
    	for col = 1 : P.n
	        xc = x[col]
	        @inbounds  for k = P.colptr[col] : (P.colptr[col+1]-1)
	            y[rv[k]] += xc / (1 << (-3 * nzv[k])) # = xc * 2 ^ (3 * nzv[k]); nzv[k] <= 0
        	end
		end
	else
		y = P * x
	end
	return y
end

# Interpolate a vector x from the global mesh to the local mesh
function interpGlobalToLocal(x::Vector{Float64}, P, y0::Vector{Float64})
	return y0 + interpGlobalToLocal(x, P)
end

# Interpolate a vector x from the global mesh to the local mesh
function interpGlobalToLocal(x::Vector{Float64}, P::SparseMatrixCSC)
	n = div(length(x), size(P,1))
	(length(x) == n * size(P,1)) || error("Incompatible input sizes")
	x = reshape(x,(size(P,1),n))
	if (eltype(P.nzval) == Int16) || (eltype(P.nzval) == Int8)
		nzv = P.nzval
		rv  = P.rowval
		y   = zeros(P.n,n)
		@inbounds begin
			for k = 1:n
				for i = 1 : P.n
					tmp = 0.0
					for j = P.colptr[i] : (P.colptr[i+1]-1)
						tmp += x[rv[j],k] / (1 << (-3 * nzv[j])) # = x[rv[j],k] * 2 ^ (3 * nzv[j]); nzv[j] <= 0
					end
					y[i,k] = tmp
				end
			end
		end
	else
		y = P' * x
	end
	return vec(y)
end
