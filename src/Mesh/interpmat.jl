export interpmat, getInterpolationMatrix
export getNodalInterpolationMatrix, getEdgeInterpolationMatrix, getFaceInterpolationMatrix, getCellCenteredInterpolationMatrix

"""
function jInv.Mesh.getInterpolationMatrix

computes bi/trilinear interpolation matrix P for cell-centered data
from Mesh1 to Mesh2. If I1 is a cell-centerd discretization of some function
on Mesh1 then, its interpolant on Mesh2 is given by

I2 = P*I1

Required Input:

	M1::AbstractTensorMesh
	M2::AbstractTensorMesh

    Note: kwargs are not used for AbstractTensorMesh

Example:

	In mesh-decoupling, we use different meshes for the inverse solution
	and the different forward problems.

	Mesh2Mesh = getInterpolationMatrix(Minv,Mfor)

"""
getInterpolationMatrix(M1::AbstractTensorMesh, M2::AbstractTensorMesh; kwargs...) =
							getCellCenteredInterpolationMatrix(M1, getCellCenteredGrid(M2))

function interpmat(x,y,Xr)
	nx = length(x); ny = length(y); np = size(Xr,1)

	# get cell sizes
	h1 = x[2:end] - x[1:end-1]
	h2 = y[2:end] - y[1:end-1]

	# find valid points
	valid1 = (x[1]-h1[1] .<=Xr[:,1]) .& (Xr[:,1] .<= x[end]+h1[end])
	valid2 = (y[1]-h2[1] .<=Xr[:,2]) .& (Xr[:,2] .<= y[end]+h2[end])
	valid =  findall( valid1 .& valid2  )

	if isempty(valid)
	    return spzeros(np,nx*ny)
	end

	J1,V1 = inter1D(Xr[:,1],x,valid)
	J2,V2 = inter1D(Xr[:,2],y,valid)

	IJ  = [  valid J1[:,1] J2[:,1]
	         valid J1[:,2] J2[:,1]
	         valid J1[:,1] J2[:,2]
	         valid J1[:,2] J2[:,2]
	        ]
	V = [
			V1[:,1].*V2[:,1]
			V1[:,2].*V2[:,1]
			V1[:,1].*V2[:,2]
			V1[:,2].*V2[:,2]
	]
	V  = V[LinearIndices(V)[findall(all(IJ[:,2:3] .> 0 ,dims=2))]]
	IJ = IJ[LinearIndices(IJ)[findall(all(IJ[:,2:3] .> 0 ,dims=2))],:]
	return sparse(IJ[:,1], LinearIndices((nx,ny))[CartesianIndex.(IJ[:,2],IJ[:,3])],V,np,nx*ny)
end

function interpmat(x,y,z,Xr)
	nx = length(x); ny = length(y); nz = length(z);
	np = size(Xr,1)

	# get cell sizes
	h1 = x[2:end] - x[1:end-1]
	h2 = y[2:end] - y[1:end-1]
	h3 = z[2:end] - z[1:end-1]

	# find valid points
	valid1 = (x[1]-h1[1] .<=Xr[:,1]) .& (Xr[:,1] .<= x[end]+h1[end])
	valid2 = (y[1]-h2[1] .<=Xr[:,2]) .& (Xr[:,2] .<= y[end]+h2[end])
	valid3 = (z[1]-h3[1] .<=Xr[:,3]) .& (Xr[:,3] .<= z[end]+h3[end])
	valid =  findall( valid1 .& valid2 .& valid3 )

	if isempty(valid)
	    return spzeros(np,nx*ny*nz)
	end

	J1,V1 = inter1D(Xr[:,1],x,valid)
	J2,V2 = inter1D(Xr[:,2],y,valid)
	J3,V3 = inter1D(Xr[:,3],z,valid)

	IJ  = [  valid J1[:,1] J2[:,1] J3[:,1]
	         valid J1[:,2] J2[:,1] J3[:,1]
	         valid J1[:,1] J2[:,2] J3[:,1]
	         valid J1[:,2] J2[:,2] J3[:,1]
	         valid J1[:,1] J2[:,1] J3[:,2]
	         valid J1[:,2] J2[:,1] J3[:,2]
	         valid J1[:,1] J2[:,2] J3[:,2]
	         valid J1[:,2] J2[:,2] J3[:,2]
	        ]
	V = [
			V1[:,1].*V2[:,1].*V3[:,1]
			V1[:,2].*V2[:,1].*V3[:,1]
			V1[:,1].*V2[:,2].*V3[:,1]
			V1[:,2].*V2[:,2].*V3[:,1]
			V1[:,1].*V2[:,1].*V3[:,2]
			V1[:,2].*V2[:,1].*V3[:,2]
			V1[:,1].*V2[:,2].*V3[:,2]
			V1[:,2].*V2[:,2].*V3[:,2]
	]
	V  = V[LinearIndices(V)[findall(all(IJ[:,2:4] .> 0 ,dims=2))]]
	IJ = IJ[LinearIndices(IJ)[findall(all(IJ[:,2:4] .> 0 ,dims=2))],:]
	return sparse(IJ[:,1], LinearIndices((nx,ny,nz))[CartesianIndex.(IJ[:,2],IJ[:,3],IJ[:,4])],V,np,nx*ny*nz)
end

function inter1D(x,ax,valid)
	# initialize
	J  = zeros(Int,length(valid),2)
	V  = zeros(length(valid),2)
	debit = valid
	for k=length(ax):-1:1; # go from right to left and search left neighbour
	    p =  (x[debit].-ax[k]) .>=0
		Ip = findall(p)
	    J[debit[Ip],1] .= k
		# check if right neighbour is in domain
	    if k < length(ax)
	        hk = ax[k+1]-ax[k]
	        J[debit[Ip],2] .= k+1
	        V[debit[Ip],1] .= (ax[k+1].-x[debit[Ip]])./hk
	        V[debit[Ip],2] .= (x[debit[Ip]].-ax[k])./hk
	    else
	        hk = ax[k]-ax[k-1]
	        V[debit[Ip],1] = (ax[k]+hk.-x[debit[Ip]])./hk
	    end
		debit = debit[findall(.!p)]
	end
	if !isempty(debit)
	    # check for points between ax[1]-h[1] and ax[1]
	    hk = ax[2]-ax[1]
	    P  = x[debit].-(ax[1]-hk)
	    p  = findall( P .>=0 )
	    J[debit[p],2] .= 1
	    V[debit[p],2] .= (x[debit[p],1].-(ax[1]-hk))/hk
		debit = debit[findall(P.<0)]
	end
	if !isempty(debit); error("something is wrong here!"); end;
	return J,V
end

function getNodalInterpolationMatrix(Mesh::AbstractTensorMesh, x)

	# axes
	xn,yn,zn = getNodalAxes(Mesh)

	# interpolate
	Q = interpmat(xn, yn, zn, x)

end

function getEdgeInterpolationMatrix(Mesh::AbstractTensorMesh, x)

	# axes
	xn,yn,zn = getNodalAxes(Mesh)
	xc,yc,zc = getCellCenteredAxes(Mesh)

	# interpolate
	Qx = interpmat(xc, yn, zn, x)
	Qy = interpmat(xn, yc, zn, x)
	Qz = interpmat(xn, yn, zc, x)

	Q  = blockdiag(blockdiag(Qx,Qy),Qz)

	return Q

end

function getFaceInterpolationMatrix(Mesh::AbstractTensorMesh, x)

	# axes
	xn,yn,zn = getNodalAxes(Mesh)
	xc,yc,zc = getCellCenteredAxes(Mesh)

	# interpolate
	Qx = interpmat(xn, yc, zc, x)
	Qy = interpmat(xc, yn, zc, x)
	Qz = interpmat(xc, yc, zn, x)

	Q  = blockdiag(blockdiag(Qx,Qy),Qz)

	return Q

end

function getCellCenteredInterpolationMatrix(Mesh::AbstractTensorMesh, x)
	axes = getCellCenteredAxes(Mesh)
	if Mesh.dim==2
		return interpmat(axes[1],axes[2],x)
	elseif Mesh.dim==3
		return interpmat(axes[1],axes[2],axes[3],x)
	end
end
