export getFaceAverageMatrix, getEdgeAverageMatrix, getNodalAverageMatrix
export ndgrid, meshgrid
export getDivergenceMatrix, getNodalGradientMatrix, getCurlMatrix, getNodalLaplacianMatrix,getNodalDivSigGradMatrix
export getEdgeMassMatrix, getdEdgeMassMatrix, getFaceMassMatrix,
       getdFaceMassMatrix,getNodalMassMatrix,getdNodalMassMatrix
export avN2C, avN2C_Nearest



## The matrix below properly averages Nodes to cells using 0.5,0.5
## However, we usually use this matrix the other way around: cells to nodes.
## Then, we need to decide what to do with the boundary. This matrix averages with 0.
## This matrix is used as default for all averging operators
function avN2C(n)
	# A = avN2C(n), 1D average operator
	I,J,V = SparseArrays.spdiagm_internal(0 => fill(.5,n), 1 => fill(.5,n))
	return sparse(I, J, V, n, n+1)
end


## The matrix below, uses a nearest neighbor interpolation on the boundaries, so that avN2C_Nearest'*ones = ones.
## This is important for geometric multigrid to work.
function avN2C_Nearest(n)
	# A = av(n), 1D average operator
	I,J,V = SparseArrays.spdiagm_internal(0 => fill(.5,n), 1 => fill(.5,n))
	av = sparse(I, J, V, n, n+1);
	av[1,1] = 1.0;
	av[end,end] = 1.0;
	return av;
end


function ddx(n)
# D = ddx(n), 1D derivative operator
	I,J,V = SparseArrays.spdiagm_internal(0 => fill(-1.0,n), 1 => fill(1.0,n)) 
	return sparse(I, J, V, n, n+1)
end

# --- linear operator
"""
 	function jInv.Mesh.getFaceAverageMatrix

	Returns Face-to-CellCenter average matrix from Mesh.Af.
	Matrix is constructed if Mesh.Af is empty.

	For 2D Mesh: Af = [A1 A2]
	For 3D Mesh: Af = [A1 A2 A3]

	Input:
		Mesh::Abstract Mesh

"""
function getFaceAverageMatrix(Mesh::AbstractTensorMesh; saveMat::Bool = true, avN2C::Function = avN2C)
	if isempty(Mesh.Af)
		if Mesh.dim==2
			A1 = kron(sparse(1.0I, Mesh.n[2], Mesh.n[2]),avN2C(Mesh.n[1]))
			A2 = kron(avN2C(Mesh.n[2]),sparse(1.0I, Mesh.n[1], Mesh.n[1]))
			Af = [A1 A2]
		elseif Mesh.dim==3
			I1 = sparse(1.0I, Mesh.n[1], Mesh.n[1]); 
			I2 = sparse(1.0I, Mesh.n[2], Mesh.n[2]);
			I3 = sparse(1.0I, Mesh.n[3], Mesh.n[3]);
			A1 = kron(I3,kron(I2,avN2C(Mesh.n[1])))
			A2 = kron(I3,kron(avN2C(Mesh.n[2]),I1))
			A3 = kron(avN2C(Mesh.n[3]),kron(I2,I1))
			Af = [A1 A2 A3]
		end
		if saveMat
			Mesh.Af = Af;
		end
		return Af;
	else
		return Mesh.Af
	end
end


"""
 	function jInv.Mesh.getEdgeAverageMatrix

	Returns Edge-to-CellCenter average matrix from Mesh.Ae.
	Matrix is constructed if Mesh.Ae is empty.

	For 3D Mesh: Ae = [A1 A2 A3]

	Input:
		Mesh::Abstract Mesh

"""
function getEdgeAverageMatrix(Mesh::AbstractTensorMesh; saveMat::Bool = true, avN2C::Function = avN2C)
	if isempty(Mesh.Ae)
		if Mesh.dim==3
			Av1 = avN2C(Mesh.n[1]);
			Av2 = avN2C(Mesh.n[2]);
			Av3 = avN2C(Mesh.n[3]);
			A1 = kron(Av3,kron(Av2,speye(Mesh.n[1]) ))
			A2 = kron(Av3,kron(speye(Mesh.n[2]),Av1))
			A3 = kron(speye(Mesh.n[3]),kron(Av2,Av1))
			Ae = [A1 A2 A3]
		elseif Mesh.dim==2
			Av1 = avN2C(Mesh.n[1]);
			Av2 = avN2C(Mesh.n[2]);
			A1 = kron(Av2,sparse(1.0I, Mesh.n[1], Mesh.n[1]))
			A2 = kron(sparse(1.0I, Mesh.n[2], Mesh.n[2]),Av1)
			Ae = [A1 A2]
		else
			error("getEdgeAverageMatrix not implemented fot $(Mesh.dim)D Meshes")
		end
		if saveMat
			Mesh.Ae = Ae;
		end
		return Ae;
	else
		return Mesh.Ae
	end
end


"""
 	function jInv.Mesh.getNodalAverageMatrix

	Returns Nodal-to-CellCenter average matrix from Mesh.An.
	Matrix is constructed if Mesh.An is empty.

	Input:
		Mesh::Abstract Mesh

"""
function getNodalAverageMatrix(Mesh::AbstractTensorMesh; saveMat::Bool = true, avN2C::Function = avN2C)
# Mesh.An = getNodalAverageMatrix(Mesh::TensorMesh3D) builds nodal-to-cc average operator
	if isempty(Mesh.An)
		if Mesh.dim==2
			An = kron(avN2C(Mesh.n[2]),avN2C(Mesh.n[1]))
		elseif Mesh.dim==3
			An = kron(avN2C(Mesh.n[3]),kron(avN2C(Mesh.n[2]),avN2C(Mesh.n[1])))
		else
			error("getNodalAverageMatrix: Dimension must be 2 or 3")
		end
		if saveMat
			Mesh.An = An;
		end
		return An;
	else
		return Mesh.An
	end
end

"""
	function jInv.Mesh.getEdgeMassMatrix(mesh,sigma)

	Returns mass matrix on cell edges, weighted by vector sigma.
	Matrix is always constructed. Uses pre-constructed edge averaging
	and cell volume matrices if available.

	Input:
		Mesh::Abstract Mesh
	       sigma::Vector

	Output:
		SparseMatrixCSC{Float64,Int64}
"""
function getEdgeMassMatrix(M::AbstractMesh,sigma::Vector)

	Ae   = getEdgeAverageMatrix(M)
	V    = getVolume(M)
	Masse    = sdiag(Ae'*(V*vec(sigma)))
	return Masse
end

"""
	function jInv.Mesh.getdEdgeMassMatrix(M,v)

	Returns directional derivative of edge mass matrix w.r.t. sigma, i.e.,

			d_sigma (M(sigma)*v)

	Matrix is always constructed. Uses pre-constructed edge averaging
	and cell volume matrices if available.

	Input:
		Mesh::Abstract - Mesh
		   v::Vector   - edge vector defining derivative

	Output:
		SparseMatrixCSC{Float64,Int64}
"""
function getdEdgeMassMatrix(M::AbstractMesh,v::Vector)

	Ae   = getEdgeAverageMatrix(M)
	V    = getVolume(M)
	return sdiag(v)*Ae'*V
end

# for compatibility with applications supporting isotropic/anisotropic PDE coefficients
getdEdgeMassMatrix(M::AbstractMesh,sigma::Vector,v::Vector) = getdEdgeMassMatrix(M, v)

"""
	function jInv.Mesh.getFaceMassMatrix(M,sigma)

	Returns face mass matrix, weighted by vector sigma.
	Matrix is always constructed. Uses pre-constructed face averaging
	and cell volume matrices if available.

	Input:
		M::Abstract Mesh
	       sigma::Vector


	Output:
		SparseMatrixCSC{Float64,Int64}
"""
function getFaceMassMatrix(M::AbstractMesh,sigma::Vector; saveMat::Bool = true, avN2C::Function = avN2C)
	Af    = getFaceAverageMatrix(M;saveMat=saveMat,avN2C=avN2C)
	V     = getVolume(M;saveMat=saveMat)
	Massf = sdiag(Af'*(V*sigma));
	return Massf
end

"""
	function jInv.Mesh.getdFaceMassMatrix(M,v)

	Returns directional derivative of face mass matrix w.r.t sigma, i.e.

		d_sigma (M(sigma)*v)

	Matrix is always constructed. Uses pre-constructed face averaging
	and cell volume matrices if available.

	Input:
		M::Abstract  - Mesh
		   v::Vector -  face vector defining derivative

	Output:
		SparseMatrixCSC{Float64,Int64}
"""
function getdFaceMassMatrix(M::AbstractMesh,v::Vector)
	Af    = getFaceAverageMatrix(M)
	V     = getVolume(M)
	Massf = sdiag(v)*Af'*V
	return Massf
end

# for compatibility with applications supporting isotropic/anisotropic PDE coefficients
getdFaceMassMatrix(M::AbstractMesh,sigma::Vector,v::Vector) = getdFaceMassMatrix(M, v)

"""
	function jInv.Mesh.getNodalMassMatrix(M,sigma)

	Returns nodal mass matrix, weighted by vector sigma.
	Matrix is always constructed. Uses pre-constructed nodal averaging
	and cell volume matrices if available.

	Input:
		M::Abstract Mesh
	       sigma::Vector

	Output:
		SparseMatrixCSC{Float64,Int64}
"""
function getNodalMassMatrix(M::AbstractMesh,sigma::Vector)
    An = getNodalAverageMatrix(M)
    V  = getVolume(M)
    Mn = Diagonal(An'*(V*sigma))
  return Mn
end

"""
	function jInv.Mesh.getdNodalMassMatrix(M,v)

	Returns directional derivative of nodal mass matrix w.r.t. sigma, i.e.,

		d_sigma (M(sigma)*v)

	Matrix is always constructed. Uses pre-constructed nodal averaging
	and cell volume matrices if available.

	Input:
		M::Abstract  - Mesh
		   v::Vector - nodal vector defining derivative

	Output:
		SparseMatrixCSC{Float64,Int64}
"""
function getdNodalMassMatrix(M::AbstractMesh,v::Vector)
	An  = getNodalAverageMatrix(M)
	V   = getVolume(M)
	dMn = Diagonal(v)*An'*V
end

# for compatibility with applications supporting isotropic/anisotropic PDE coefficients
getdNodalMassMatrix(M::AbstractMesh,sigma::Vector,v::Vector) = getdNodalMassMatrix(M, v)

# --- ndgrid
ndgrid(v::AbstractVector) = copy(v)
# function ndgrid{T}(v1::AbstractVector{T}, v2::AbstractVector{T})
# 	m, n = length(v1), length(v2)
# 	v1 = reshape(v1, m, 1)
# 	v2 = reshape(v2, 1, n)
# 	(repeat(v1, 1, n), repeat(v2, m, 1))
# end

function ndgrid_fill(a, v, s, snext)
	for j = 1:length(a)
		a[j] = v[div(rem(j-1, snext), s)+1]
	end
end

function ndgrid(vs::AbstractVector{T}...) where {T}
	n = length(vs)
	sz = map(length, vs)
	out = ntuple(i->Array{T}(undef,sz), n)
	s = 1
	for i=1:n
		a = out[i]::Array
		v = vs[i]
		snext = s*size(a,i)
		ndgrid_fill(a, v, s, snext)
		s = snext
	end
	out
end

# --- meshgrid
meshgrid(v::AbstractVector) = meshgrid(v, v)
function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T}) where {T}
	m, n = length(vy), length(vx)
	vx = reshape(vx, 1, n)
	vy = reshape(vy, m, 1)
	(repeat(vx, m, 1), repeat(vy, 1, n))
end

function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T},
	                 vz::AbstractVector{T}) where {T}
	m, n, o = length(vy), length(vx), length(vz)
	vx = reshape(vx, 1, n, 1)
	vy = reshape(vy, m, 1, 1)
	vz = reshape(vz, 1, 1, o)
	om = ones(Int, m)
	on = ones(Int, n)
	oo = ones(Int, o)
	(vx[om, :, oo], vy[:, on, oo], vz[om, on, :])
end

# --- differential operators
function getDivergenceMatrix(Mesh::AbstractTensorMesh;saveMat = true)
# Mesh.Div = getDivergenceMatrix(Mesh::AbstractTensorMesh) builds face-to-cc divergence operator
	if isempty(Mesh.Div)
		if Mesh.dim==2
			D1 = kron(sparse(1.0I, Mesh.n[2], Mesh.n[2]),ddx(Mesh.n[1]))
			D2 = kron(ddx(Mesh.n[2]),sparse(1.0I, Mesh.n[1], Mesh.n[1]))
			Div = [D1 D2]
		elseif Mesh.dim==3
			D1 = kron(sparse(1.0I, Mesh.n[3], Mesh.n[3]),kron(sparse(1.0I, Mesh.n[2], Mesh.n[2]),ddx(Mesh.n[1])))
			D2 = kron(sparse(1.0I, Mesh.n[3], Mesh.n[3]),kron(ddx(Mesh.n[2]),sparse(1.0I, Mesh.n[1], Mesh.n[1])))
			D3 = kron(ddx(Mesh.n[3]),kron(sparse(1.0I, Mesh.n[2], Mesh.n[2]),sparse(1.0I, Mesh.n[1], Mesh.n[1])))
			Div = [D1 D2 D3]
		end
		Vi = getVolumeInv(Mesh)
		F  = getFaceArea(Mesh)
		Div = Vi*(Div*F);
		if saveMat
			Mesh.Div = Div;
		end
		return Div;
	else
		return Mesh.Div;
	end
end

function getNodalGradientMatrix(Mesh::AbstractTensorMesh; saveMat::Bool = true)
# Mesh.Grad = getNodalGradientMatrix(Mesh::AbstractTensorMesh) builds nodal-to-edge gradient operator
	if isempty(Mesh.Grad) 
		if Mesh.dim==2
			G1 = kron(sparse(1.0I, Mesh.n[2]+1, Mesh.n[2]+1),ddx(Mesh.n[1]))
			G2 = kron(ddx(Mesh.n[2]),sparse(1.0I, Mesh.n[1]+1, Mesh.n[1]+1))
			Grad =[G1; G2]
		elseif Mesh.dim==3
			G1 = kron(sparse(1.0I, Mesh.n[3]+1, Mesh.n[3]+1),kron(sparse(1.0I, Mesh.n[2]+1, Mesh.n[2]+1),ddx(Mesh.n[1])))
			G2 = kron(sparse(1.0I, Mesh.n[3]+1, Mesh.n[3]+1),kron(ddx(Mesh.n[2]),sparse(1.0I, Mesh.n[1]+1, Mesh.n[1]+1)))
			G3 = kron(ddx(Mesh.n[3]),kron(sparse(1.0I, Mesh.n[2]+1, Mesh.n[2]+1),sparse(1.0I, Mesh.n[1]+1, Mesh.n[1]+1)))
			Grad =[G1; G2; G3]
		end
		Li  = getLengthInv(Mesh)
		Grad = Li*Grad;
		if saveMat
			Mesh.Grad = Grad
		end
		return Grad;
	else
		return Mesh.Grad
	end
	
end

function getCurlMatrix(Mesh::AbstractTensorMesh)
# Mesh.Curl = getCurlMatrix(Mesh::AbstractTensorMesh) builds edge-to-face curl operator
	if isempty(Mesh.Curl)
		if Mesh.dim==3
			# The Curl from edges to faces
			Dyz = kron(ddx(Mesh.n[3]),kron(sparse(1.0I, Mesh.n[2], Mesh.n[2]),sparse(1.0I, Mesh.n[1]+1, Mesh.n[1]+1)))
			Dzy = kron(sparse(1.0I, Mesh.n[3], Mesh.n[3]),kron(ddx(Mesh.n[2]),sparse(1.0I, Mesh.n[1]+1, Mesh.n[1]+1)))

			Dxz = kron(ddx(Mesh.n[3]),kron(sparse(1.0I, Mesh.n[2]+1, Mesh.n[2]+1),sparse(1.0I, Mesh.n[1], Mesh.n[1])))
			Dzx = kron(sparse(1.0I, Mesh.n[3], Mesh.n[3]),kron(sparse(1.0I, Mesh.n[2]+1, Mesh.n[2]+1),ddx(Mesh.n[1])))

			Dxy = kron(sparse(1.0I, Mesh.n[3]+1, Mesh.n[3]+1),kron(ddx(Mesh.n[2]),sparse(1.0I, Mesh.n[1], Mesh.n[1])))
			Dyx = kron(sparse(1.0I, Mesh.n[3]+1, Mesh.n[3]+1),kron(sparse(1.0I, Mesh.n[2], Mesh.n[2]),ddx(Mesh.n[1])))

			# curl on the edges
			Curl = [
				 spzeros(Mesh.nf[1],Mesh.ne[1])  -Dyz   Dzy
				 Dxz   spzeros(Mesh.nf[2],Mesh.ne[2])  -Dzx
				-Dxy   Dyx   spzeros(Mesh.nf[3],Mesh.ne[3])]

				Fi = getFaceAreaInv(Mesh)
				L  = getLength(Mesh)
				Mesh.Curl = Fi*(Curl*L)
		else
			error("CURL is only implemented for 3D meshes.")
		end

	end
	return Mesh.Curl
end

function getNodalLaplacianMatrix(Mesh::AbstractTensorMesh; saveMat::Bool = true)
	
	if isempty(Mesh.nLap)
		G = getNodalGradientMatrix(Mesh)
		L = G'*G;
		if saveMat 
			Mesh.nLap = L;
		end
		return L;
	else
		return Mesh.nLap;
	end
	
	
end

function getNodalDivSigGradMatrix(M::AbstractTensorMesh,sig::Vector; saveMat::Bool = true,avN2C::Function = avN2C)
 	G       = getNodalGradientMatrix(M, saveMat = saveMat)
	Ae      = getEdgeAverageMatrix(M, saveMat = saveMat,avN2C = avN2C);
    L       = (G'*sparse(Diagonal(Ae'*(sig))))*G
	return L
end
