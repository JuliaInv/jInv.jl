export getEdgeIntegralOfPolygonalChain,getStraightLineCurrentIntegral

#using Base.BLAS

"""
        function jInv.Mesh.getEdgeIntegralOfPolygonalChain

        s = getEdgeIntegralPolygonalChain(mesh,polygon)
        s = getEdgeIntegralPolygonalChain(mesh,polygon,normalize)

        Compute the integral of a piecewise linear edge grid basis projected onto
        the edges of a polygonal chain. This function can be used to evaluate the
        source term of a line current carried by the polygonal chain in electromagnetic
        modelling.

        The piecewise linear edge grid basis consists of the functions
          phix_i (i = 1 ... numXEdges)
          phiy_i (i = 1 ... numYEdges)
          phiz_i (i = 1 ... numZEdges)
        where phix_i, phiy_i, phiz_i = 1 at the i-th x/y/z-directed edge
          and phix_i, phiy_i, phiz_i = 0 at all other edges.

        INPUT
        mesh ...... Tensor mesh
        polygon ... vertices of polygonal chain as numVertices x 3 array
        normalize . divide line integral by length of integration path (boolean, optional)

        OUTPUT
        s ......... source vector

        For a closed current loop, specify polygon such that
        polygon[1,:] == polygon[end,:]
"""
function getEdgeIntegralOfPolygonalChain(M::TensorMesh3D,polygon::Array{Float64,2};
                                         normalize=false)

#Reformat input
nx, ny, nz =  M.n[1],  M.n[2],  M.n[3]
x0, y0, z0 = M.x0[1], M.x0[2], M.x0[3]
hx, hy, hz = M.h1, M.h2, M.h3
px, py, pz = polygon[:,1], polygon[:,2], polygon[:,3]

# nodal grid
x,y,z = getNodalAxes(M)
# = [x0; x0 + cumsum(hx)]
# y  = [y0; y0 + cumsum(hy)]
# z  = [z0; z0 + cumsum(hz)]

# discrete edge function
sx = zeros(nx  , ny+1, nz+1)
sy = zeros(nx+1, ny  , nz+1)
sz = zeros(nx+1, ny+1, nz  )

# number of line segments
np = length(px) - 1

# check that all polygon vertices are inside the mesh
for ip = 1:np+1
  ax = px[ip]
  ay = py[ip]
  az = pz[ip]
  ix = findlast((ax .>= x[1:nx-1]) .& (ax .<= x[2:nx]))
  iy = findlast((ay .>= y[1:ny-1]) .& (ay .<= y[2:ny]))
  iz = findlast((az .>= z[1:nz-1]) .& (az .<= z[2:nz]))
  if (ix < 1) | (iy < 1) | (iz < 1)
    msg = @sprintf("Polygon vertex (%d,%d,%d) is outside the mesh \n", ax, ay, az)
    error(msg)
  end
end

# integrate each line segment
for ip = 1:np

    # start and end vertices
    ax = px[ip]
    ay = py[ip]
    az = pz[ip]
    bx = px[ip+1]
    by = py[ip+1]
    bz = pz[ip+1]

    # find intersection with mesh planes
    dx = bx - ax
    dy = by - ay
    dz = bz - az
    d  = sqrt(dx^2 + dy^2 + dz^2)
    tol = d * eps(Float64)
    if abs(dx) > tol
      tx = (x .- ax) / dx
      tx = tx[(tx .>=0) .& (tx .<= 1)]
    else
      tx = Float64[]
    end
    if abs(dy) > tol
      ty = (y .- ay) / dy
      ty = ty[(ty .>=0) .& (ty .<= 1)]
    else
      ty = Float64[]
    end
    if abs(dz) > tol
      tz = (z .- az) / dz
      tz = tz[(tz .>=0) .& (tz .<= 1)]
    else
      tz = Float64[]
    end

    t  = sort(unique([-0.0;0.0;tx;ty;tz;1.0]))[2:end]
    nq = length(t) - 1
    tc = 0.5 * (t[1:nq] + t[2:nq+1])

    for iq = 1:nq

        cx = ax + tc[iq] * dx
        cy = ay + tc[iq] * dy
        cz = az + tc[iq] * dz

        # locate cell id
        ix = findlast((cx .>= x[1:nx-1]) .& (cx .<= x[2:nx]))
        iy = findlast((cy .>= y[1:ny-1]) .& (cy .<= y[2:ny]))
        iz = findlast((cz .>= z[1:nz-1]) .& (cz .<= z[2:nz]))

        # local coordinates
        hxloc = hx[ix]
        hyloc = hy[iy]
        hzloc = hz[iz]
        axloc = ax + t[iq]   * dx - x[ix]
        ayloc = ay + t[iq]   * dy - y[iy]
        azloc = az + t[iq]   * dz - z[iz]
        bxloc = ax + t[iq+1] * dx - x[ix]
        byloc = ay + t[iq+1] * dy - y[iy]
        bzloc = az + t[iq+1] * dz - z[iz]
        # integrate
        sxloc,syloc,szloc = getStraightLineCurrentIntegral(hxloc,hyloc,hzloc,axloc,ayloc,azloc,bxloc,byloc,bzloc)

        # integrate
        sx[ix,[iy,iy+1],[iz,iz+1]] = sx[ix,[iy,iy+1],[iz,iz+1]] + reshape(sxloc, (2, 2))
        sy[[ix,ix+1],iy,[iz,iz+1]] = sy[[ix,ix+1],iy,[iz,iz+1]] + reshape(syloc, (2, 2))
        sz[[ix,ix+1],[iy,iy+1],iz] = sz[[ix,ix+1],[iy,iy+1],iz] + reshape(szloc, (2, 2))
    end
end

s = [sx[:]; sy[:]; sz[:]]

# normalize
if normalize
	if all(polygon[1,:] .== polygon[np+1,:])

		# closed polygon: divide by enclosed area
		a  = 0.0
		px = polygon[2:np,1] .- polygon[1,1]
		py = polygon[2:np,2] .- polygon[1,2]
		pz = polygon[2:np,3] .- polygon[1,3]
		for ip = 1:np-2
			cx = py[ip] * pz[ip+1] - pz[ip] * py[ip+1]
			cy = pz[ip] * px[ip+1] - px[ip] * pz[ip+1]
			cz = px[ip] * py[ip+1] - py[ip] * px[ip+1]
			a += sqrt(cx * cx + cy * cy + cz * cz)
		end
		a *= 0.5

	else
		# open polygon: divide by length
		a = 0.0
		for ip = 1:np
			dx = polygon[ip+1,1] - polygon[ip,1]
			dy = polygon[ip+1,2] - polygon[ip,2]
			dz = polygon[ip+1,3] - polygon[ip,3]
			a += sqrt(dx * dx + dy * dy + dz * dz)
		end

	end

	s /= a

end

return s

end

"""
        function jInv.Mesh.getStraightLineCurrentIntegral

        [sx,sy,sz] = getStraightLineCurrentIntegral(hx,hy,hz,ax,ay,az,bx,by,bz)

        Compute integral int(W . J dx^3) in brick of size hx x hy x hz
        where W denotes the 12 local bilinear edge basis functions
        and where J prescribes a unit line current
        between points (ax,ay,az) and (bx,by,bz).
"""
function getStraightLineCurrentIntegral(hx,hy,hz,ax,ay,az,bx,by,bz)

  # length of line segment
  lx = bx - ax
  ly = by - ay
  lz = bz - az
  l  = sqrt(lx^2 + ly^2 + lz^2)

  if l == 0
    sx = zeros(4,1)
    sy = zeros(4,1)
    sz = zeros(4,1)
    return
  end

  # linear interpolation between a and b
  x(t)=ax + t * lx  #0 <= t <= 1
  y(t)=ay + t * ly
  z(t)=az + t * lz

  # edge basis functions
  wx(t)=[
    (1. - y(t) / hy) * (1. - z(t) / hz)
    (    y(t) / hy) * (1. - z(t) / hz)
    (1. - y(t) / hy) * (    z(t) / hz)
    (    y(t) / hy) * (    z(t) / hz)]
  wy(t) = [
    (1. - x(t) / hx) * (1. - z(t) / hz)
    (    x(t) / hx) * (1. - z(t) / hz)
    (1. - x(t) / hx) * (    z(t) / hz)
    (    x(t) / hx) * (    z(t) / hz)]
  wz(t) = [
    (1. - x(t) / hx) * (1. - y(t) / hy)
    (    x(t) / hx) * (1. - y(t) / hy)
    (1. - x(t) / hx) * (    y(t) / hy)
    (    x(t) / hx) * (    y(t) / hy)]

  # integration using Simpson's rule
  sx = (wx(0.) + 4. * wx(0.5) + wx(1.)) * (lx / 6.)
  sy = (wy(0.) + 4. * wy(0.5) + wy(1.)) * (ly / 6.)
  sz = (wz(0.) + 4. * wz(0.5) + wz(1.)) * (lz / 6.)

  return sx, sy, sz
end
