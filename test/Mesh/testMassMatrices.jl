using Test
using jInv.Mesh
using jInv.Utils


n      = [8 15 13]

# regular mesh
domain = [0.1 1.1 -1.2 1.4 0.5 1.3]
Mr     = getRegularMesh(domain,n)

# tensor mesh
h1     = rand(n[1])
h2     = rand(n[2])
h3     = rand(n[3])
Mt     = getTensorMesh3D(h1,h2,h3)

# conductivity
sig    = rand(prod(n))
Meshes = (Mr,Mt)

massMatrices = ( (getEdgeMassMatrix,getdEdgeMassMatrix),      
                 (getNodalMassMatrix,getdNodalMassMatrix), 
                 (getFaceMassMatrix,getdFaceMassMatrix))

for k=1:length(Meshes)
	for j=1:length(massMatrices)
		println("\ttest $(massMatrices[j][1]) on $(typeof(Meshes[k]))")
		M  = Meshes[k]
		Mf3 = massMatrices[j][1]
		dMf3 = massMatrices[j][2]
		Av = Mf3(M,sig);
		vt = randn(size(Av,1))
		
	
		chkDir, = checkDerivative(sig->Mf3(M,sig)*vt,sig->dMf3(M,sig,vt),sig,out=false)
		@test chkDir
	end
end
