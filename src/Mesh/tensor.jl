export TensorMesh3D, getTensorMesh3D
export getCellCenteredGrid, getNodalGrid, getFaceGrids, getEdgeGrids
export getCellCenteredAxes, getNodalAxes,getBoundaryNodes
export getVolume, getVolumeInv, getFaceArea, getFaceAreaInv, getLength, getLengthInv

"""
	mutable struct jInv.Mesh.TensorMesh3D <: AbstractTensorMesh
		
	Fields:
	
		h1::Vector{Float64}     - cell size in x1 direction
		h2::Vector{Float64}     - cell size in x2 direction
		h3::Vector{Float64}     - cell size in x3 direction
		x0::Vector{Float64}     - origin
		dim::Int                - dimension (dim=3)
		n::Vector{Int64}        - number of cells in each direction
		nc::Int                 - nc total number of cells (nc=prod(n))
		nf::Vector{Int64}       - number of faces
		ne::Vector{Int64}       - number of edges
	
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
	
	Example: 
	
	h1 = rand(4); h2 = rand(6); h3 = rand(5);
	M  = getTensorMesh2D(h1,h2,h3)
	
"""
mutable struct TensorMesh3D <: AbstractTensorMesh
	h1::Vector{Float64}
	h2::Vector{Float64}
	h3::Vector{Float64}
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
	function jInv.Mesh.getTensorMesh3D
		
	constructs TensorMesh3D
	
	Required Input:
	
	   h1::Array - cell-sizes in x1 direction
	   h2::Array - cell-sizes in x2 direction
	   h3::Array - cell-sizes in x3 direction

	Optional Input
	   x0::Array - origin (default = zeros(3))

"""
function getTensorMesh3D(h1::Array,h2::Array,h3::Array,x0=zeros(3))
	n = [length(h1); length(h2);  length(h3)]
	nc = prod(n)
	nf = [(n[1]+1)*n[2]*n[3]; n[1]*(n[2]+1)*n[3]; n[1]*n[2]*(n[3]+1) ]
	ne = [n[1]*(n[2]+1)*(n[3]+1); (n[1]+1)*n[2]*(n[3]+1); (n[1]+1)*(n[2]+1)*n[3] ]
	empt = spzeros(0,0);
	dim = 3
return  TensorMesh3D(h1,h2,h3,x0,dim,n,nc,nf,ne,empt,empt,empt,empt,empt,empt,empt,empt,empt,empt,empt,empt,empt)
end

import Base.==
function ==(M1::TensorMesh3D,M2::TensorMesh3D)
	isEqual = fill(true,20)

	# check mandatory fields
	isEqual[1] =  (M1.h1 == M2.h1) & (M1.h2 == M2.h2) & (M1.h3 == M2.h3)
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


function getNodalAxes(Mesh::TensorMesh3D)
	xn,yn,zn = getNodalAxes(Mesh.h1,Mesh.h2,Mesh.h3)
	return xn.+Mesh.x0[1],yn.+Mesh.x0[2],zn.+Mesh.x0[3]
end

function getCellCenteredAxes(Mesh::TensorMesh3D)
	xc,yc,zc = getCellCenteredAxes(Mesh.h1,Mesh.h2,Mesh.h3)
	return xc.+Mesh.x0[1],yc.+Mesh.x0[2],zc.+Mesh.x0[3]
end

function getNodalGrid(Mesh::TensorMesh3D)
# X = getNodalGrid(Mesh::TensorMesh3D)
	xn,yn,zn = getNodalAxes(Mesh)
	Xn,Yn,Zn = ndgrid(xn,yn,zn)
	return [vec(Xn) vec(Yn) vec(Zn)]
end

function getCellCenteredGrid(Mesh::TensorMesh3D)
# X = getCellCenteredGrid(Mesh::TensorMesh3D)
	xc,yc,zc = getCellCenteredAxes(Mesh)
	Xc,Yc,Zc = ndgrid(xc,yc,zc);
	return [vec(Xc) vec(Yc) vec(Zc)]
end

function getEdgeGrids(Mesh::TensorMesh3D)
# [Xe1, Xe2, Xe3] = getEdgeGrids(Mesh::TensorMesh3D)

	xn,yn,zn    = getNodalAxes(Mesh)
	xc,yc,zc    = getCellCenteredAxes(Mesh)
	
	Xe1,Ye1,Ze1 = ndgrid(xc,yn,zn)
	Xe2,Ye2,Ze2 = ndgrid(xn,yc,zn)
	Xe3,Ye3,Ze3 = ndgrid(xn,yn,zc)
	
	return  [vec(Xe1) vec(Ye1) vec(Ze1)], [vec(Xe2) vec(Ye2) vec(Ze2)], [vec(Xe3) vec(Ye3) vec(Ze3)]

end

function getFaceGrids(Mesh::TensorMesh3D)
# [Xf1, Xf2, Xf3] = getFaceGrids(Mesh::TensorMesh3D)

	xn,yn,zn    = getNodalAxes(Mesh)
	xc,yc,zc    = getCellCenteredAxes(Mesh)
	
	Xf1,Yf1,Zf1 = ndgrid(xn,yc,zc)
	Xf2,Yf2,Zf2 = ndgrid(xc,yn,zc)
	Xf3,Yf3,Zf3 = ndgrid(xc,yc,zn)
	
	return [vec(Xf1) vec(Yf1) vec(Zf1)], [ vec(Xf2) vec(Yf2) vec(Zf2)], [vec(Xf3) vec(Yf3) vec(Zf3)]	

end


function getNodalAxes(h1,h2,h3)
	nc = [length(h1); length(h2);  length(h3)]

	x = zeros(nc[1]+1); for i=1:nc[1], x[i+1] = x[i] + h1[i]; end
	y = zeros(nc[2]+1); for i=1:nc[2], y[i+1] = y[i] + h2[i]; end
	z = zeros(nc[3]+1); for i=1:nc[3], z[i+1] = z[i] + h3[i]; end
	return x,y,z
end

function getCellCenteredAxes(h1,h2,h3)
	x,y,z = getNodalAxes(h1,h2,h3)
	# cell centered grids
	xc = x[1:end-1] + h1/2;
	yc = y[1:end-1] + h2/2;
	zc = z[1:end-1] + h3/2;
	return xc,yc,zc
end

function getNodalGrid(h1,h2,h3)
# X = getNodalGrid(h1,h2,h3)
	xn,yn,zn = getNodalAxes(h1,h2,h3)
	Xn,Yn,Zn = ndgrid(xn,yn,zn)
	return [vec(Xn) vec(Yn) vec(Zn)]

end

function getCellCenteredGrid(h1,h2,h3)
# X = getCellCenteredGrid(h1,h2,h3)
	xc,yc,zc = getCellCenteredAxes(h1,h2,h3)
	Xc,Yc,Zc = ndgrid(xc,yc,zc);
	return [vec(Xc) vec(Yc) vec(Zc)]
end

function getEdgeGrids(h1,h2,h3)
# [Xe1,Xe2,Xe3] = getEdgeGrids(h1,h2,h3)
	xn,yn,zn = getNodalAxes(h1,h2,h3)
	xc,yc,zc = getCellCenteredAxes(h1,h2,h3)
	
	Xe1,Ye1,Ze1 = ndgrid(xc,yn,zn)
	Xe2,Ye2,Ze2 = ndgrid(xn,yc,zn)
	Xe3,Ye3,Ze3 = ndgrid(xn,yn,zc)
	return  [vec(Xe1) vec(Ye1) vec(Ze1)], [vec(Xe2) vec(Ye2) vec(Ze2)], [vec(Xe3) vec(Ye3) vec(Ze3)]

end

function getFaceGrids(h1,h2,h3)
# [Xf1, Xf2, Xf3] = getFaceGrids(h1,h2,h3)
	xn,yn,zn = getNodalAxes(h1,h2,h3)
	xc,yc,zc = getCellCenteredAxes(h1,h2,h3)
	
	Xf1,Yf1,Zf1 = ndgrid(xn,yc,zc)
	Xf2,Yf2,Zf2 = ndgrid(xc,yn,zc)
	Xf3,Yf3,Zf3 = ndgrid(xc,yc,zn)
	
	return [vec(Xf1) vec(Yf1) vec(Zf1)], [ vec(Xf2) vec(Yf2) vec(Zf2)], [vec(Xf3) vec(Yf3) vec(Zf3)]	
end

# --- linear operators for tensor mesh
function getVolume(Mesh::TensorMesh3D;saveMat::Bool=true)
# Mesh.V = getVolume(Mesh::TensorMesh3D) computes volumes v, returns diag(v)
	if isempty(Mesh.Vi)
		V = kron(sdiag(Mesh.h3),kron(sdiag(Mesh.h2),sdiag(Mesh.h1)))
		if saveMat
			Mesh.V = V;
		end
		return V;
	else
		return Mesh.V
	end
end
function getVolumeInv(Mesh::TensorMesh3D;saveMat::Bool=true)
# Mesh.Vi = getVolumeInv(Mesh::TensorMesh3D) returns sdiag(1 ./v)
	if isempty(Mesh.Vi)
		Vi = kron(sdiag(1 ./Mesh.h3),kron(sdiag(1 ./Mesh.h2),sdiag(1 ./Mesh.h1)))
		if saveMat
			Mesh.Vi = Vi;
		end
		return Vi;
	else
		return Mesh.Vi
	end
end
function getFaceArea(Mesh::TensorMesh3D)
# Mesh.F = getFaceArea(Mesh::TensorMesh3D) computes face areas a, returns  sdiag(a)
	if isempty(Mesh.F)
		f1  = kron(sdiag(Mesh.h3)      ,kron(sdiag(Mesh.h2)      ,sparse(1.0I,Mesh.n[1]+1,Mesh.n[1]+1)))
		f2  = kron(sdiag(Mesh.h3)      ,kron(sparse(1.0I,Mesh.n[2]+1,Mesh.n[2]+1) ,sdiag(Mesh.h1)))
		f3  = kron(sparse(1.0I,Mesh.n[3]+1,Mesh.n[3]+1) ,kron(sdiag(Mesh.h2)      ,sdiag(Mesh.h1)))
		Mesh.F = blockdiag(blockdiag(f1,f2),f3)
	end
return Mesh.F
end
function getFaceAreaInv(Mesh::TensorMesh3D)
# Mesh.Fi = getFaceAreaInv(Mesh::TensorMesh3D) computes inverse of face areas, returns sdiag(1 ./a)
	if isempty(Mesh.Fi)
		f1i  = kron(sdiag(1 ./Mesh.h3)   ,kron(sdiag(1 ./Mesh.h2)   ,sparse(1.0I,Mesh.n[1]+1,Mesh.n[1]+1)))
		f2i  = kron(sdiag(1 ./Mesh.h3)   ,kron(sparse(1.0I,Mesh.n[2]+1,Mesh.n[2]+1) ,sdiag(1 ./Mesh.h1)))
		f3i  = kron(sparse(1.0I,Mesh.n[3]+1,Mesh.n[3]+1) ,kron(sdiag(1 ./Mesh.h2)   ,sdiag(1 ./Mesh.h1)))
		Mesh.Fi = blockdiag(blockdiag(f1i,f2i),f3i)
	end
return Mesh.Fi
end

function getLength(Mesh::TensorMesh3D)
# Mesh.L = getLength(Mesh::TensorMesh3D) computes edge lengths l, returns sdiag(l)
	if isempty(Mesh.L)
		l1  = kron(sparse(1.0I,Mesh.n[3]+1,Mesh.n[3]+1),kron(sparse(1.0I,Mesh.n[2]+1,Mesh.n[2]+1),sdiag(Mesh.h1)))
		l2  = kron(sparse(1.0I,Mesh.n[3]+1,Mesh.n[3]+1),kron(sdiag(Mesh.h2),sparse(1.0I,Mesh.n[1]+1,Mesh.n[1]+1)))
		l3  = kron(sdiag(Mesh.h3),kron(sparse(1.0I,Mesh.n[2]+1,Mesh.n[2]+1),sparse(1.0I,Mesh.n[1]+1,Mesh.n[1]+1)))
		Mesh.L   = blockdiag(blockdiag(l1,l2),l3);
	end
return Mesh.L
end

function getLengthInv(Mesh::TensorMesh3D)
# Mesh.L = getLength(Mesh::TensorMesh3D) computes inverse of edge lengths l, returns sdiag(1 ./l)
	if isempty(Mesh.Li)
		l1i = kron(sparse(1.0I,Mesh.n[3]+1,Mesh.n[3]+1),kron(sparse(1.0I,Mesh.n[2]+1,Mesh.n[2]+1),sdiag(1 ./Mesh.h1)))
		l2i = kron(sparse(1.0I,Mesh.n[3]+1,Mesh.n[3]+1),kron(sdiag(1 ./Mesh.h2),sparse(1.0I,Mesh.n[1]+1,Mesh.n[1]+1)))
		l3i = kron(sdiag(1 ./Mesh.h3),kron(sparse(1.0I,Mesh.n[2]+1,Mesh.n[2]+1),sparse(1.0I,Mesh.n[1]+1,Mesh.n[1]+1)))
		Mesh.Li  = blockdiag(blockdiag(l1i,l2i),l3i);
	end
return Mesh.Li
end

"""
function jInv.Mesh.getBoundaryNodes(Mesh)
	
	Returns an array tuple (boundary node indices, inner node indices)
	
	Input: 
		Mesh::Abstract Mesh
	
	Output:
		Tuple{Array{Int64,1},Array{Int64,1}}
"""
function getBoundaryNodes(M::AbstractTensorMesh)
    
    # number of cells in each direction
    N = M.n
    
    # nodal matrix
    nodal_mat = 1:prod(N.+1)
    
    # check dimension
    if M.dim==2
        
        nodal_mat = reshape(nodal_mat, N[1]+1, N[2]+1)
        
        iin = nodal_mat[2:end-1,2:end-1]
        ib = setdiff(nodal_mat, iin)
        
    elseif M.dim==3
        
        nodal_mat = reshape(nodal_mat, N[1]+1, N[2]+1, N[3]+1)
        
        iin = nodal_mat[2:end-1,2:end-1,2:end-1]
        ib = setdiff(nodal_mat, iin)
        
    end
    
    return (vec(ib), vec(iin))
end

include("getEdgeIntegralOfPolygonalChain.jl")
