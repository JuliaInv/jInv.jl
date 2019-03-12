using jInv.InverseSolve

n = 100 * 3

dc   = rand(n).+2
dobs = rand(n).+3
wd   = rand(n).+1
dd   = rand(n)

for misfitFun in [ SSDFun HuberFun ]

	# Gradient test for ", misfitFun, " (real)")

	function testFun(x,v=[])
		mis,dmis,d2mis = misfitFun(x,dobs,wd)
		if isempty(v)
			return mis
		else
			return mis,dot(dmis,v)
		end
	end

	chkDer, = checkDerivative(testFun,dc,out=false)
	@test chkDer

end

for misfitFun in [ SSDFun ]

	# Hessian test for ", misfitFun, " (real)")

	function testFun(x,v=[])
		mis,dmis,d2mis = misfitFun(x,dobs,wd)
		if isempty(v)
			return dmis
		else
			return dmis,d2mis.*v
		end
	end

	chkDer, = checkDerivative(testFun,dc,out=false)
	@test chkDer
end

### Complex

dc   += 1im*rand(n)
dobs += 1im*rand(n)
wd   += 1im*rand(n)
dd   += 1im*rand(n)

for misfitFun in [ SSDFun ]

	# Gradient test for ", misfitFun, " (complex)")

	function testFun(x,v=[])
		mis,dmis,d2mis = misfitFun(x,dobs,wd)
		if isempty(v)
			return mis
		else
			return mis,real(dot(dmis,v))
		end
	end

	chkDer, = checkDerivative(testFun,dc,out=false)
	@test chkDer

end

# for misfitFun in [ SSDFun ]

	# Hessian test for ", misfitFun, " (complex)")
	# warn("skipped"); break
# 	function testFun(x,v=[])
# 		mis,dmis,d2mis = misfitFun(x,dobs,wd)
# 		if isempty(v)
# 			return real(dmis)
# 		else
# 			return real(dmis),real(d2mis.*v)
# 		end
# 	end
#
# 	chkDer, = checkDerivative(testFun,dc,out=true)
# 	@test chkDer
#
# end
