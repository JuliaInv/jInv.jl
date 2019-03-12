using jInv.InverseSolve
using Test
using jInv.Mesh
using jInv.Utils

# build regular mesh and Iact
domain = [0;5;0;5; 0;3]
n     = [12;13;10]
M     = getRegularMesh(domain,n)
idx   = ones(Int,(n[1],n[2],n[3]))
idx[:,:,1:5] .= 0
Iact  = sparse(I,M.nc,M.nc)
Iact  = Iact[:,vec(idx).==1]
mc    = randn(size(Iact,2))

regFuns = [
			(m,mref,M)->diffusionReg(m,mref,M,Iact=Iact),
			(m,mref,M)->wdiffusionReg(m,mref,M,Iact=Iact),
			smallnessReg,
			(m,mref,M)->wTVReg(m,mref,M,Iact=Iact),
			(m,mref,M)->anisoTVReg(m,mref,M,Iact=Iact)]
for k=1:length(regFuns)
	# checkDerivative of $(regFuns[k])

	function testFun(x,v=[])
		Sc,dS,d2S = regFuns[k](x,0.0.*x,M)
		if isempty(v)
			return Sc
		else
			return Sc,dot(dS,v)
		end
	end
	chkDer, = checkDerivative(testFun,mc,out=false)
	@test chkDer
end

regFuns = [wdiffusionRegNodal;wTVRegNodal]
idx   = ones(Int,(n[1]+1,n[2]+1,n[3]+1))
idx[:,:,1:5] .= 0
Iact  = sparse(I,prod(M.n.+1),prod(M.n.+1))
Iact  = Iact[:,vec(idx).==1]
mc    = randn(size(Iact,2))

for k=1:length(regFuns)
	# checkDerivative of $(regFuns[k])

	function testFun(x,v=[])
		Sc,dS,d2S = regFuns[k](x,0.0.*x,M,Iact=Iact)
		if isempty(v)
			return Sc
		else
			return Sc,dot(dS,v)
		end
	end
	chkDer, = checkDerivative(testFun,mc,out=false)
	@test chkDer
end


regFuns = [logBarrier;logBarrierSquared]
idx   = ones(Int,(n[1]+1,n[2]+1,n[3]+1))
idx[:,:,1:5] .= 0
Iact  = sparse(I,prod(M.n.+1),prod(M.n.+1))
Iact  = Iact[:,vec(idx).==1]
mc    = ones(size(Iact,2)) + 0.05*randn(size(Iact,2));

low   = ones(length(mc))*0.4;
high  = ones(length(mc))*1.6;
epsilon = ones(length(mc))*0.58;

for k=1:length(regFuns)
	# checkDerivative of $(regFuns[k])

	function testFun(x,v=[])
		Sc,dS,d2S = regFuns[k](x,0.0.*x,M,low,high,epsilon)
		if isempty(v)
			return Sc
		else
			return Sc,dot(dS,v)
		end
	end
	chkDer, = checkDerivative(testFun,mc,out=false)
	@test chkDer
end
