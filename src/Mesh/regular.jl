export RegularMesh, getRegularMesh
export getCellCenteredGrid, getNodalGrid, getEdgeGrids, getFaceGrids
export getCellCenteredAxes, getNodalAxes
export getVolume, getVolumeInv, getFaceArea, getFaceAreaInv, getLength, getLengthInv

"""
	mutable struct jInv.Mesh.RegularMesh <: AbstractTensorMesh
		
	Regular mesh in 1D, 2D, and 3D
	
	Fields:
		
		domain::Vector{Float64}  - physical domain [min(x1) max(x1) min(x2) max(x2)]
		h::Vector{Float64}       - cell size
		x0::Vector{Float64}      - origin
		dim::Int                 - dimension of mesh
		n::Vector{Int64}         - number of cells in each dimension
		nc::Int                  - total number of cells
		nf::Vector{Int64}        - number of faces in each dimension
		ne::Vector{Int64}        - number of edges in each dimension
		
		
		Persistent Operators:

		Operators should not be accessed directly. They will be built, if needed,
		when accessing them using specified method. clear!(M) will release all 
		memory.

			Div::SparseMatrixCSC    - divergence (faces -> cell-centers)
									  Access via: getDivergenceMatrix(M)
			Grad::SparseMatrixCSC   - gradient (nodal -> edges)
									  Access via: getNodalGradientMatrix(M)
			Curl::SparseMatrixCSC   - curl (edges -> faces)
									  Access via: getCurlMatrix(M)
			Af::SparseMatrixCSC     - face average (faces -> cell-centers)
									  Access via: getFaceAverageMatrix(M)
			Ae::SparseMatrixCSC     - edge average (edges -> cell-centers)
									  Access via: getEdgeAverageMatrix(M)
			An::SparseMatrixCSC     - nodal average (nodes -> cell-centers)
									  Access via: getNodalAverageMatrix(M)
			V::SparseMatrixCSC      - cell volumes (diagonal matrix)
									  Access via: getVolume(M)
			F::SparseMatrixCSC      - face area (diagonal matrix)
									  Access via: getFaceArea(M)
			L::SparseMatrixCSC      - edge length (diagonal matrix)
									  Access via: getLength(M)
			Vi::SparseMatrixCSC     - inverse cell volumes (diagonal matrix)
									  Access via: getVolumeInv(M)
			Fi::SparseMatrixCSC     - inverse face area (diagonal matrix)
									  Access via: getFaceAreaInv(M)
			Li::SparseMatrixCSC     - inverse edge length (diagonal matrix)
									  Access via: getLengthAreaInv(M)
			nLap::SparseMatrixCSC   - nodal Laplacian
									  Access via: getNodalLaplacian(M)
		
	Examples:
	M2D  = getRegularMesh([1.2 2.4 2.2 5.0],[3,4])	
	M3D  = getRegularMesh([1.2 2.4 2.2 5.0 0 1],[3,4,7])	
		
"""
mutable struct RegularMesh <: AbstractTensorMesh
	domain::Vector{Float64} 
	h::Vector{Float64} 
	x0::Vector{Float64}
	dim::Int
	n::Vector{Int64}
	nc::Int
	nf::Vector{Int64}
	ne::Vector{Int64}
	Div::SparseMatrixCSC
	Grad::SparseMatrixCSC
	Curl::SparseMatrixCSC
	Af::SparseMatrixCSC
	Ae::SparseMatrixCSC
	An::SparseMatrixCSC
	V::SparseMatrixCSC
	F::SparseMatrixCSC
	L::SparseMatrixCSC
	Vi::SparseMatrixCSC
	Fi::SparseMatrixCSC
	Li::SparseMatrixCSC
	nLap::SparseMatrixCSC
end

"""
	function jInv.Mesh.getRegularMesh
		
	Constructs regular mesh
	
	Input: 
		domain - physical domain rectangular
		n      - number of cells in each dimension

	Examples:
	M2D  = getRegularMesh([1.2 2.4 2.2 5.0],[3,4])	
	M3D  = getRegularMesh([1.2 2.4 2.2 5.0 0 1],[3,4,7])	

"""
function getRegularMesh(domain,n)
	domain = vec(float(domain))
	n     = vec(n)
	nc    = prod(n)
	h     = vec((domain[2:2:end]-domain[1:2:end])./n)
	dim   = round(Int64,length(domain)/2)
	if dim==1
		nf = [(n[1]+1)]
		ne = [(n[1]+1)]
	elseif dim==2
		nf = [(n[1]+1)*n[2]; n[1]*(n[2]+1);]
		ne = [n[1]*(n[2]+1); (n[1]+1)*n[2];  ]
	elseif dim==3
		nf = [(n[1]+1)*n[2]*n[3]; n[1]*(n[2]+1)*n[3]; n[1]*n[2]*(n[3]+1); ]
		ne = [n[1]*(n[2]+1)*(n[3]+1); (n[1]+1)*n[2]*(n[3]+1); (n[1]+1)*(n[2]+1)*n[3]; ]
	end
	x0 = domain[1:2:end]
	empt = spzeros(0,0);
return RegularMesh(domain,h,x0,dim,n,nc,nf,ne,empt,empt,empt,empt,empt,empt,empt,empt,empt,empt,empt,empt,empt)
end


import Base.==
function ==(M1::RegularMesh,M2::RegularMesh)
	isEqual = fill(true,20)
	
	# check mandatory fields
	isEqual[1] =  (M1.domain == M2.domain)
	isEqual[2] =  (M1.h     == M2.h)
	isEqual[3] =  (M1.x0    == M2.x0)
	isEqual[4] =  (M1.dim   == M2.dim)
	isEqual[5] =  (M1.n     == M2.n)
	isEqual[6] =  (M1.nc    == M2.nc)
	isEqual[7] =  (M1.nf    == M2.nf)
	isEqual[8] =  (M1.ne    == M2.ne)
	
	# check fields that might be empty
	if !(isempty(M1.Div)) && !(isempty(M2.Div))
		isEqual[9] = (M1.Div == M2.Div)
	end
	if !(isempty(M1.Grad)) && !(isempty(M2.Grad))
		isEqual[10] = (M1.Grad == M2.Grad)
	end
	if !(isempty(M1.Curl)) && !(isempty(M2.Curl))
		isEqual[11] = (M1.Curl == M2.Curl)
	end
	if !(isempty(M1.Af)) && !(isempty(M2.Af))
		isEqual[12] = (M1.Af == M2.Af)
	end
	if !(isempty(M1.Ae)) && !(isempty(M2.Ae))
		isEqual[13] = (M1.Ae == M2.Ae)
	end
	if !(isempty(M1.An)) && !(isempty(M2.An))
		isEqual[14] = (M1.An == M2.An)
	end
	if !(isempty(M1.V)) && !(isempty(M2.V))
		isEqual[15] = (M1.V == M2.V)
	end
	if !(isempty(M1.F)) && !(isempty(M2.F))
		isEqual[16] = (M1.F == M2.F)
	end
	if !(isempty(M1.L)) && !(isempty(M2.L))
		isEqual[17] = (M1.L == M2.L)
	end
	if !(isempty(M1.Vi)) && !(isempty(M2.Vi))
		isEqual[18] = (M1.Vi == M2.Vi)
	end
	if !(isempty(M1.Fi)) && !(isempty(M2.Fi))
		isEqual[19] = (M1.Fi == M2.Fi)
	end
	if !(isempty(M1.Li)) && !(isempty(M2.Li))
		isEqual[20] = (M1.Li == M2.Li)
	end
	return all(isEqual)
end


# --- grid constructor
function getCellCenteredGrid(Mesh::RegularMesh)
# X = getCellCenteredGrid(Mesh::RegularMesh)
	return getCellCenteredGrid(Mesh.domain,Mesh.n)
end

function getNodalGrid(Mesh::RegularMesh)
# X = getNodalGrid(Mesh::RegularMesh)
	return getNodalGrid(Mesh.domain,Mesh.n)
end

function getEdgeGrids(Mesh::RegularMesh)
	return getEdgeGrids(Mesh.domain,Mesh.n)
end

function getFaceGrids(Mesh::RegularMesh)
	return getFaceGrids(Mesh.domain,Mesh.n)
end

function getCellCenteredGrid(domain,n)
# X = getCellCenteredGrid(domain,n)
	dim = round(Int64,length(domain)/2)
	if dim==1
		xc = getCellCenteredAxes(domain,n)
	elseif dim==2
		x1,x2 = getCellCenteredAxes(domain,n)
		X1,X2 = ndgrid(x1,x2)
		xc = [vec(X1) vec(X2)]
	elseif dim==3
		x1,x2,x3 = getCellCenteredAxes(domain,n)
		X1,X2,X3 = ndgrid(x1,x2,x3)
		xc = [vec(X1) vec(X2) vec(X3)]
	end
	return xc
end

function getNodalGrid(domain,n)
# X = getNodalGrid(domain,nc)
	dim = round(Int64,length(domain)/2)
	if dim==1
		xc = getNodalAxes(domain,n)
	elseif dim==2
		x1,x2 = getNodalAxes(domain,n)
		X1,X2 = ndgrid(x1,x2)
		xc = [vec(X1) vec(X2)]
	elseif dim==3
		x1,x2,x3 = getNodalAxes(domain,n)
		X1,X2,X3 = ndgrid(x1,x2,x3)
		xc = [vec(X1) vec(X2) vec(X3)]
	end
	return xc
end

function getEdgeGrids(domain,nc)
# X = getEdgeGrids(domain,nc)
	dim = round(Int64,length(domain)/2)
	h   = (domain[2:2:end]-domain[1:2:end])./nc
	if dim==2
		x1n,x2n = getNodalAxes(domain,nc)
		x1c,x2c = getCellCenteredAxes(domain,nc)
		# edge-1 grid
		X1t,X2t = ndgrid(x1c,x2n)
		x1 = [vec(X1t) vec(X2t)]
		# edge-2 grid
		X1t,X2t = ndgrid(x1n,x2c)
		x2 = [vec(X1t) vec(X2t)]
		return (x1,x2)
	elseif dim==3
		x1n,x2n,x3n = getNodalAxes(domain,nc)
		x1c,x2c,x3c = getCellCenteredAxes(domain,nc)
		# edge-1 grid
		X1,X2,X3 = ndgrid(x1c,x2n,x3n)
		x1 = [vec(X1) vec(X2) vec(X3)]
		# edge-2 grid
		X1,X2,X3 = ndgrid(x1n,x2c,x3n)
		x2 = [vec(X1) vec(X2) vec(X3)]
		# edge-3 grid
		X1,X2,X3 = ndgrid(x1n,x2n,x3c)
		x3 = [vec(X1) vec(X2) vec(X3)]
		return (x1,x2,x3)
	end
end

function getFaceGrids(domain,nc)
# X = getFaceGrids(domain,nc)
	dim = round(Int64,length(domain)/2)
	h   = (domain[2:2:end]-domain[1:2:end])./nc
	if dim==2
		x1n,x2n = getNodalAxes(domain,nc)
		x1c,x2c = getCellCenteredAxes(domain,nc)
		# face-1 grid
		X1t,X2t = ndgrid(x1n,x2c)
		x1 = [vec(X1t) vec(X2t)]
		# face-2 grid
		X1t,X2t = ndgrid(x1c,x2n)
		x2 = [vec(X1t) vec(X2t)]
		return (x1,x2)
	elseif dim==3
		x1n,x2n,x3n = getNodalAxes(domain,nc)
		x1c,x2c,x3c = getCellCenteredAxes(domain,nc)
		# face-1 grid
		X1,X2,X3 = ndgrid(x1n,x2c,x3c)
		x1 = [vec(X1) vec(X2) vec(X3)]
		# face-2 grid
		X1,X2,X3 = ndgrid(x1c,x2n,x3c)
		x2 = [vec(X1) vec(X2) vec(X3)]
		# face-3 grid
		X1,X2,X3 = ndgrid(x1c,x2c,x3n)
		x3 = [vec(X1) vec(X2) vec(X3)]
		return (x1,x2,x3)
	end
end

function getNodalAxes(Mesh::RegularMesh)
	return getNodalAxes(Mesh.domain,Mesh.n)
end

function getNodalAxes(domain,nc)
	dim = round(Int64,length(domain)/2)
	
	if dim==1
		x1 = collect(range(domain[1], step=(domain[2]-domain[1])/nc,length=nc+1))
		return x1
	elseif dim==2
		x1 = collect(range(domain[1], step=(domain[2]-domain[1])/nc[1],length=nc[1]+1))
		x2 = collect(range(domain[3], step=(domain[4]-domain[3])/nc[2],length=nc[2]+1))
		return (x1,x2)
	elseif dim==3
		x1 = collect(range(domain[1], step=(domain[2]-domain[1])/nc[1],length=nc[1]+1))
		x2 = collect(range(domain[3], step=(domain[4]-domain[3])/nc[2],length=nc[2]+1))
		x3 = collect(range(domain[5], step=(domain[6]-domain[5])/nc[3],length=nc[3]+1))
		return (x1,x2,x3)
	end
end

function getCellCenteredAxes(Mesh::RegularMesh)
	return getCellCenteredAxes(Mesh.domain,Mesh.n)
end

function getCellCenteredAxes(domain,nc)
	dim = round(Int64,length(domain)/2)
	h   = vec(domain[2:2:end]-domain[1:2:end])./vec(nc)
	if dim==1
		x1 = collect(range(domain[1]+h[1]/2, step=h[1], length=nc[1]))
		return x1
	elseif dim==2
		x1 = collect(range(domain[1]+h[1]/2, step=h[1], length=nc[1]))
		x2 = collect(range(domain[3]+h[2]/2, step=h[2], length=nc[2]))
		return (x1,x2)
	elseif dim==3
		x1 = collect(range(domain[1]+h[1]/2, step=h[1], length=nc[1]))
		x2 = collect(range(domain[3]+h[2]/2, step=h[2], length=nc[2]))
		x3 = collect(range(domain[5]+h[3]/2, step=h[3], length=nc[3]))
		return (x1,x2,x3)
	end
end

# --- linear operators for tensor mesh
function getVolume(Mesh::RegularMesh;saveMat::Bool=true)
# Mesh.V = getVolume(Mesh::RegularMesh) computes volumes v, returns diag(v)
	if isempty(Mesh.V)
		V = sparse(prod(Mesh.h)*I, Mesh.nc,Mesh.nc)
		if saveMat
			Mesh.V = V;
		end
		return V;
	else
		return Mesh.V;
	end
end

function getVolumeInv(Mesh::RegularMesh; saveMat::Bool = true)
# Mesh.Vi = getVolumeInv(Mesh::RegularMesh) returns sdiag(1 ./v)
	if isempty(Mesh.Vi)
		Vi = sparse((1/prod(Mesh.h))*I, Mesh.nc,Mesh.nc)
		if saveMat
			Mesh.Vi = Vi;
		end
		return Vi;
	else
		return Mesh.Vi;
	end
end

function getFaceArea(Mesh::RegularMesh; saveMat = true)
# Mesh.F = getFaceArea(Mesh::RegularMesh) computes face areas a, returns  sdiag(a)
	if isempty(Mesh.F)
		if Mesh.dim==2
			f1  = Mesh.h[2]*sparse(1.0I,Mesh.nf[1],Mesh.nf[1])
			f2  = Mesh.h[1]*sparse(1.0I,Mesh.nf[2],Mesh.nf[2])
			F = blockdiag(f1,f2)
		elseif Mesh.dim==3
			f1  = (Mesh.h[3]*Mesh.h[2])*sparse(1.0I,Mesh.nf[1],Mesh.nf[1])
			f2  = (Mesh.h[3]*Mesh.h[1])*sparse(1.0I,Mesh.nf[2],Mesh.nf[2])
			f3  = (Mesh.h[2]*Mesh.h[1])*sparse(1.0I,Mesh.nf[3],Mesh.nf[3])
			F = blockdiag(blockdiag(f1,f2),f3)
		end
		if saveMat
			Mesh.F = F;
		end
		return F;
	else
		return Mesh.F;
	end
end
function getFaceAreaInv(Mesh::RegularMesh)
# Mesh.Fi = getFaceAreaInv(Mesh::RegularMesh) computes inverse of face areas, returns sdiag(1 ./a)
	if isempty(Mesh.Fi)
		if Mesh.dim==2
			f1i  = (1/Mesh.h[2])*sparse(1.0I,Mesh.nf[1],Mesh.nf[1])
			f2i  = (1/Mesh.h[1])*sparse(1.0I,Mesh.nf[2],Mesh.nf[2])
			Mesh.Fi = blockdiag(f1i,f2i)
		elseif Mesh.dim==3
			f1i  = (1/(Mesh.h[3]*Mesh.h[2]))*sparse(1.0I,Mesh.nf[1],Mesh.nf[1])
			f2i  = (1/(Mesh.h[3]*Mesh.h[1]))*sparse(1.0I,Mesh.nf[2],Mesh.nf[2])
			f3i  = (1/(Mesh.h[2]*Mesh.h[1]))*sparse(1.0I,Mesh.nf[3],Mesh.nf[3])
			Mesh.Fi = blockdiag(blockdiag(f1i,f2i),f3i)
		end
	end
	return Mesh.Fi
end

function getLength(Mesh::RegularMesh)
# Mesh.L = getLength(Mesh::RegularMesh) computes edge lengths l, returns sdiag(l)
	if isempty(Mesh.L)
		if Mesh.dim==2
			l1  = Mesh.h[1]*sparse(1.0I,Mesh.ne[1],Mesh.ne[1])
			l2  = Mesh.h[2]*sparse(1.0I,Mesh.ne[2],Mesh.ne[2])
			Mesh.L   = blockdiag(l1,l2)
		elseif Mesh.dim==3
			l1  = Mesh.h[1]*sparse(1.0I,Mesh.ne[1],Mesh.ne[1])
			l2  = Mesh.h[2]*sparse(1.0I,Mesh.ne[2],Mesh.ne[2])
			l3  = Mesh.h[3]*sparse(1.0I,Mesh.ne[3],Mesh.ne[3])
			Mesh.L   = blockdiag(blockdiag(l1,l2),l3)
		end
	end
	return Mesh.L
end

function getLengthInv(Mesh::RegularMesh)
# Mesh.L = getLength(Mesh::RegularMesh) computes inverse of edge lengths l, returns sdiag(1 ./l)
	if isempty(Mesh.Li)
		if Mesh.dim==2
			l1i  = (1/Mesh.h[1])*sparse(1.0I,Mesh.ne[1],Mesh.ne[1])
			l2i  = (1/Mesh.h[2])*sparse(1.0I,Mesh.ne[2],Mesh.ne[2])
			Mesh.Li   = blockdiag(l1i,l2i)
		elseif Mesh.dim==3
			l1i  = (1/Mesh.h[1])*sparse(1.0I,Mesh.ne[1],Mesh.ne[1])
			l2i  = (1/Mesh.h[2])*sparse(1.0I,Mesh.ne[2],Mesh.ne[2])
			l3i  = (1/Mesh.h[3])*sparse(1.0I,Mesh.ne[3],Mesh.ne[3])
			Mesh.Li  = blockdiag(blockdiag(l1i,l2i),l3i)
		end
	end
	return Mesh.Li
end


