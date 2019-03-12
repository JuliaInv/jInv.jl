using Distributed
using SparseArrays
using Printf

import jInv.InverseSolve
import jInv.Mesh
import jInv.LinearSolvers
import jInv.ForwardShare
import jInv.Utils

# @everywhere begin
using jInv.InverseSolve
using jInv.Mesh
using jInv.LinearSolvers
using jInv.ForwardShare
using jInv.Utils
using Test

mutable struct LSparam <: ForwardProbType
	A
	Ainv
end

import jInv.ForwardShare.getNumberOfData
function getNumberOfData(pFor::LSparam)
	return size(pFor.A,1)
end
import jInv.ForwardShare.getSensMatSize
function getSensMatSize(pFor::LSparam)
	return size(pFor.A)
end

import jInv.ForwardShare.getData
function getData(m::Vector,pFor::LSparam,doClear::Bool=false)
	  d = pFor.A*m
	  if doClear
	    clear!(pFor)
          end
	  return d,pFor
end

import jInv.ForwardShare.getSensMatVec
function getSensMatVec(v::Vector,m::Vector,pFor::LSparam)
	return pFor.A*v
end

import jInv.ForwardShare.getSensTMatVec
function getSensTMatVec(v::Vector,m::Vector,pFor::LSparam)
	return pFor.A'*v
end

import jInv.Utils.clear!
function clear!(pFor::LSparam)
	pFor.A = SparseMatrixCSC(I, 0, 0);
	pFor.Ainv = [];
end

# end
