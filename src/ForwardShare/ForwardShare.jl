"""
jInv.ForwardShare

This module contains generic methods for forward problem computations. It handles the
parallelisms and abstractions used in jInv to allow fast and flexible solution of
PDE parameter estimation problems.

Most users will not use ForwardShare directly, but extend its methods by customized
implementations of their problems. Doing so, users can benefit from the full parallel
functionality.

Here is a quick step-by-step introuction on how to define your own forward problem in
a way that is compatible with jInv.

1) Create your own type (documented param that describes a forward problem) as a subtype
   of the generic ForwardProbType

   		type YourParam <: ForwardProbType
      		Mesh
      		Sources
      		Receivers
          Fields
	  	end

2) Import and extend the method 'getData' for your type. Implement it such that it solves
   one forward problem that is stored locally on a worker. The I/O should be as follows

	    function jInv.ForwardShare.getData(sigma,param::YourParam)
        # solve your problem
	    	return dobs, param
	    end

3) Import and extend methods for computing matrix vector products with the sensitivity
   matrix for your problem.

      function jInv.ForwardShare.getSensMatVec(x,sigma,param::YourParam)
        # compute stuff
   			return Sens*x
      end

      function jInv.ForwardShare.getSensTMatVec(x,sigma,param::YourParam)
        # compute stuff
      	return transpose(Sens) * x
      end

	4) clear! method for your type to get rid of unnecessary temp results (e.g., fields, factorizations, preconditioners )

  5) With the previous 4 steps you are already in good shape using ForwardShare to
     reconstruct parameters in parallel or in combination with other physical models.
     For some advanced functionality you will need to extend the methods

        jInv.ForwardShare.getNumberOfData
        jInv.ForwardShare.getSensMatSize
"""
module ForwardShare


	export ForwardProbType
	abstract type ForwardProbType end

	using jInv.Mesh
	using jInv.Utils
	
	using Distributed
	using SparseArrays
	using LinearAlgebra

	export getSensMatVec
	"""
	Jv  = getSensMatVec(v::Vector,m::Vector,param::ForwardProbType)

	Computes matrix-vector product with the Jacobian.

	"""
	getSensMatVec(v::Vector,m::Vector,param::ForwardProbType) = error("nyi")
	export getSensTMatVec
	"""
	JTv  = getSensMatVec(v::Vector,m::Vector,param::ForwardProbType)

	Computes matrix-vector product with the transpose of Jacobian. Implementation
	depends on forward problem.

	"""
	getSensTMatVec(v::Vector,m::Vector,param::ForwardProbType) = error("nyi")

	"""
	(m,n) = getSensMatSize(pFor)

	Returns size of sensitivity matrix where m is the number of data points
	and n the number of parameters in the model.

	Input

		pFor - forward problem:: Union{ForwardProbType, Array, RemoteChannel}

	This is problem dependent and should be implemented in the respective
	packages.
	"""
	getSensMatSize(pFor::ForwardProbType) = error("getSensMatSize not implemented for forward problems of type $(typeof(pFor))")

	function getSensMatSize(pFor::RemoteChannel)
		if pFor.where != myid()
			return remotecall_fetch(getSensMatSize,pFor.where,pFor)
		else
			return getSensMatSize(fetch(pFor))
		end
	end

	function getSensMatSize(pFor::Array)
		n  = length(pFor)
		sz = [0;0]
		for k=1:n
			(s1,s2) = getSensMatSize(pFor[k])
			sz[2]=s2
			sz[1]+=s1
		end
		return tuple(sz...)
	end

  """
  nd = getNumberOfData(pFor)

  Returns number of data in forward problem
  """
  getNumberOfData(pFor::ForwardProbType) = error("getNumberOfData not implemented for forward problems of type $(typeof(pFor))")

  function getNumberOfData(pFor::RemoteChannel)
		if pFor.where != myid()
			return remotecall_fetch(getNumberOfData,pFor.where,pFor)
		else
			return getNumberOfData(fetch(pFor))
		end
	end

	function getNumberOfData(pFor::Array)
		n  = length(pFor)
		nd = 0
		for k=1:n
			nd += getNumberOfData(pFor[k])
		end
		return nd
	end

	export getSensTMatVec,getSensMatVec, getSensMatSize, getNumberOfData

	Mesh2MeshTypes = Union{RemoteChannel, Future, SparseMatrixCSC, AbstractFloat}

	# # ===== Methods for parallelization =====
	include("getDataParallel.jl")
	include("prepareMesh2Mesh.jl")
	include("interpLocalToGlobal.jl")
	include("getSensMat.jl")
	include("testing.jl")

	import jInv.Utils.clear!
	function clear!(P::ForwardProbType;clearAinv::Bool=true,clearFields::Bool=true, clearMesh::Bool=false, clearSources::Bool=false, clearObs::Bool=false,clearAll::Bool=false)
		if clearAll || clearMesh
			P.M = clear!(P.M)
		end
		if clearAll || clearSources
			P.Sources = clear!(P.Sources)
		end
		if clearAll || clearObs
			P.Obs     = clear!(P.Obs)
		end
		if clearAll || clearFields
		if any(fieldnames(typeof(P)) .== :Fields)
			P.Fields = clear!(P.Fields)
		end

	end

	if clearAll || clearAinv
	if any(fieldnames(typeof(P)) .== :Ainv)
	   clear!(P.Ainv)
	end
	end


end


end
