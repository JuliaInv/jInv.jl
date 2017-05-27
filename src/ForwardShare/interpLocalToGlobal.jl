export interpLocalToGlobal, interpGlobalToLocal
# Interpolate a vector x from the local mesh to the global mesh

function interpLocalToGlobal(x::Vector{Float64}, P::AbstractFloat,y0::AbstractFloat=0.0)
	return P * x + y0
end

function interpGlobalToLocal(x::Vector{Float64}, P::AbstractFloat,y0::AbstractFloat=0.0)
	return P * x + y0 
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
	if (eltype(P.nzval) == Int16) || (eltype(P.nzval) == Int8)
		nzv = P.nzval
		rv  = P.rowval
		y   = zeros(P.n)  
		@inbounds begin
			for i = 1 : P.n
				tmp = 0.0
				for j = P.colptr[i] : (P.colptr[i+1]-1)
					tmp += x[rv[j]] / (1 << (-3 * nzv[j])) # = x[rv[j]] * 2 ^ (3 * nzv[j]); nzv[j] <= 0
				end
				y[i] = tmp
			end
		end
	else
		y = P' * x
	end
	return y
end