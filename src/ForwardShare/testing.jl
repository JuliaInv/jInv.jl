export adjointTest

function adjointTest(m0,pFor,out::Bool=true)
	
end
	

import jInv.Utils.checkDerivative
function checkDerivative(m0,pFor::ForwardProbType;kwargs...)
    function testFun(x,v=[])
        fc = vec(getData(x,pFor)[1])
		if !(isempty(v))
			return fc, getSensMatVec(v,x,pFor)
		else
			return fc
		end
	end
    return checkDerivative(testFun,m0;kwargs...)
end

"""
function jInv.adjointTest
	
automatic adjoint test for forward problems. If sensitivity matrix is m x n then
this method generates vector v,m of lengths n and m, respectively and computes

err = abs(dot(v,JTw) - dot(w,Jv))  / abs(dot(w,Jv))

Input:

	sig::Vector           - current model
	pFor::ForwardProbType - forward problem
	

Optional Inputs:

	out::Bool            - controls verbosity (default=false)
	tol::Real            - tol on relative error (default=1e-10)
	
Output: 

	passed              - true/false depending whether test passed or not
	err                 - absolute error
"""
function adjointTest(sig,pFor::ForwardProbType;out::Bool=false,tol::Real=1e-10)
    
    (out) && println("calling forward problem")
    dobs, pFor = getData(sig,pFor)
    
    (out) && println("calling getSensMatVec")
    v = getRandomTestDirection(sig)
    Jv = getSensMatVec(v,sig,pFor)
    w = getRandomTestDirection(Jv)
    t1 = dot(w,Jv)
    (out) && println("t1 = $(t1)")
    
    (out) && println("calling getSensTMatVec")
    JTw = getSensTMatVec(w,sig,pFor)
    t2 = dot(v,JTw)
    (out) && println("t2 = $(t2)")
    
    err = abs(t1-t2)/abs(t1)
    return err < tol, err
end
