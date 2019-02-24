using jInv.InverseSolve
using Test
using jInv.Mesh
using jInv.Utils

# setup forward mesh
domain = [0;5;0;5; 0;3]
n     = [12;13;10]
Minv  = getRegularMesh(domain,n)

modFun = (expMod,boundMod,fMod)
for k=1:length(modFun)
	# checkDerivative of $(modFun[k])
	function testModFun(m,v=[])
		mc,dm = modFun[k](m)
		if !isempty(v)
			dm = dm*v
			return mc,dm
		end
		return mc
	end
	derivativeDistModel, = checkDerivative(testModFun,rand(prod(n)),out=false)
	@test derivativeDistModel
end
