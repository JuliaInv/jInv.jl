module Utils

using Printf
using Distributed
using SparseArrays
using LinearAlgebra

	include("sparseUtils.jl")
	include("testing.jl")
	include("expandPolygon.jl")
	include("sortpermFast.jl")
	include("uniqueidx.jl")
	include("initRemoteChannel.jl")
	include("variousUtils.jl")



  import Distributed.clear!
	export clear!

	# function Base.sub2ind(n::Array{Int64,1},ii::Array{Int64,1},jj::Array{Int64,1},kk::Array{Int64,1})
	# 	return Base.sub2ind((n[1],n[2],n[3]),ii,jj,kk)
	# end
	#
	# function Base.sub2ind(n::Array{Int64,1},ii::Int64,jj::Int64,kk::Int64)
	# 	return Base.sub2ind((n[1],n[2],n[3]),ii,jj,kk)
	# end

	function clear!(R::RemoteChannel)
		p = take!(R)
		p = clear!(p)
		put!(R,p);
		return R;
	end

	function clear!(F::Future)
		p = fetch(F)
		p = clear!(p)
		put!(F,p);
		return F;
	end

	function clear!(PF::Array{RemoteChannel})
		@sync begin
			for p=workers()
				@async begin
					for i=1:length(PF)
						if p==PF[i].where
							PF[i] = initRemoteChannel(clear!, p, PF[i])
						end
					end
				end
			end
		end
	end

	function clear!(PF::Array{Future})
		@sync begin
			for p=workers()
				@async begin
					for i=1:length(PF)
						if p==PF[i].where
							PF[i] = remotecall(clear!, p, PF[i])
						end
					end
				end
			end
		end
	end


	function clear!(x::Array{T,N}) where {T,N}
		return Array{T,N}(undef,ntuple((i)->0, N))
	end

	# function clear!{T}(x::Vector{T})
		# return Array(T,0)
	# end

	function clear!(A::SparseMatrixCSC{T}) where {T}
		return spzeros(0,0);
	end

	export getWorkerIds
	function getWorkerIds(A::Array{RemoteChannel})
		Ids = []
		for k=1:length(A)
			push!(Ids,A[k].where)
		end
		return unique(Ids)
	end
end
