using jInv.Mesh
using Test

# setup 3D problem
nc = [6, 7, 8]
x0 = rand(3)
domain = [x0[1], 4, x0[2], 2, x0[3], 6]
h   = (domain[2:2:end]-domain[1:2:end])./nc
h1  = h[1]*ones(nc[1])
h2  = h[2]*ones(nc[2])
h3  = h[3]*ones(nc[3])

Mt            = getTensorMesh3D(h1,h2,h3,x0)
l             = diag(getLength(Mt))
nex           = Mt.ne[1]
ne            = sum(Mt.ne)
Xe1,Xe2,Xe3   = getEdgeGrids(h1,h2,h3)
Xe1           = Xe1 .+ x0'

xt1 = x0[1] +   h[1]
xt2 = x0[1] + 5*h[1]
yt1 = x0[2] + 3*h[2]
yt2 = x0[2] + 4*h[2]
zt1 = x0[3] + 4*h[3]
zt2 = x0[3] + 5*h[3]

# test polygonal path directly on edges
idx     = findall( (xt1 .< Xe1[:,1] .< xt2) .& (Xe1[:,2] .== yt1) .& (Xe1[:,3] .== zt1) )
s1      = zeros(ne)
s1[idx] = l[idx]
p1      = [xt1 yt1 zt1
           xt2 yt1 zt1]
s1Comp  = getEdgeIntegralOfPolygonalChain(Mt,p1)

@test norm(s1-s1Comp) < 1e-10

s1Norm = getEdgeIntegralOfPolygonalChain(Mt,p1,normalize=true)
@test abs(sum(s1Norm) - 1.0) < 1e-10

# test path offset by half a cell in the y direction
idx1     = findall( (xt1 .< Xe1[:,1] .< xt2) .& (Xe1[:,2] .== yt1) .& (Xe1[:,3] .== zt1) )
idx2     = findall( (xt1 .< Xe1[:,1] .< xt2) .& (Xe1[:,2] .== yt2) .& (Xe1[:,3] .== zt1) )
s2       = zeros(ne)
s2[idx1] = l[idx1]./2
s2[idx2] = l[idx2]./2
p2       = [xt1 (yt1+yt2)/2 zt1
            xt2 (yt1+yt2)/2 zt1]

s2Comp   = getEdgeIntegralOfPolygonalChain(Mt,p2)
@test norm(s2-s2Comp) < 1e-10

# test path offset by half a cell in the y and z directions
idx3     = findall( (xt1 .< Xe1[:,1] .< xt2) .& (Xe1[:,2] .== yt1) .& (Xe1[:,3] .== zt2) )
idx4     = findall( (xt1 .< Xe1[:,1] .< xt2) .& (Xe1[:,2] .== yt2) .& (Xe1[:,3] .== zt2) )
s3       = zeros(ne)
s3[idx1] = l[idx1]./4
s3[idx2] = l[idx2]./4
s3[idx3] = l[idx3]./4
s3[idx4] = l[idx4]./4
p3       = [xt1 (yt1+yt2)/2 (zt1 + zt2)/2
            xt2 (yt1+yt2)/2 (zt1 + zt2)/2]

s3Comp   = getEdgeIntegralOfPolygonalChain(Mt,p3)
@test norm(s3-s3Comp) < 1e-10
