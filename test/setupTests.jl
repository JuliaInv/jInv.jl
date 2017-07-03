import jInv.InverseSolve
import jInv.Mesh
import jInv.LinearSolvers
import jInv.ForwardShare
import jInv.Utils
import Base.Test
@everywhere begin
using jInv.InverseSolve
using jInv.Mesh
using jInv.LinearSolvers
using jInv.ForwardShare
using jInv.Utils
using Base.Test

mutable struct LSparam <: ForwardProbType
	A
	Ainv
end

function jInv.ForwardShare.getNumberOfData(pFor::LSparam)
	return size(pFor.A,1)
end
function jInv.ForwardShare.getSensMatSize(pFor::LSparam)
	return size(pFor.A)
end

function jInv.ForwardShare.getData(m::Vector,pFor::LSparam,doClear::Bool=false)
	  d = pFor.A*m
	  if doClear
	    clear!(pFor)
          end
	  return d,pFor
end

function jInv.ForwardShare.getSensMatVec(v::Vector,m::Vector,pFor::LSparam)
	return pFor.A*v
end

function jInv.ForwardShare.getSensTMatVec(v::Vector,m::Vector,pFor::LSparam)
	return pFor.A'*v
end

import jInv.Utils.clear!
function clear!(pFor::LSparam)
	pFor.A = speye(0);
	pFor.Ainv = [];
end

end
