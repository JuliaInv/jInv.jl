export uniqueidx, sortunique

"""
function b,ii,jj = uniqueidx(a::Vector{Int64})

b = sorted unique(a)
ii, jj are permutations:
a[ii] = b,  b[jj] = a
"""
function uniqueidx( a::Array{T,1} ) where {T<:Integer}

   n = length(a)

   #ii = sortperm(a)
   #ii = sortperm(a,alg=MergeSort)
   #b = a[ii]  # sorted array

   ii,b = sortpermFast(a)


  # jj = invperm(ii)  # inverse permutation
   jj = Array{T}(undef,n)


   nunique = n  # counter for # of unique values
   p1 = 1  # pointer to the next unique element
   jj[ii[1]] = 1

   for i = 2:n

      if b[p1] == b[i]

         iip1 = ii[p1]
         iii = ii[i]

         jj[iip1] = p1
         jj[iii]  = p1
         idx = min( iip1, iii )
         ii[p1] = idx

         nunique = nunique - 1

      else
         p1 += 1
         b[p1] = b[i]

         iii = ii[i]
         ii[p1] = iii
         jj[iii] = p1
      end

   end  # i

	deleteat!(b,  nunique+1:n)
	deleteat!(ii, nunique+1:n)

	return b, ii, jj
end # function uniqueidx

#-----------------------------------------------------

"""
function b = sortunique( a::Array{Int64,1} )

b = sorted unique(a)
same as uniqueidx but doesn't compute permutation
vectors
"""
function sortunique( a::Array{Int64,1} )

   n = length(a)

   ii,b = sortpermFast(a)
   ii = []

   nunique = n  # counter for # of unique values
   p1 = 1  # pointer to the next unique element

   for i = 2:n

      if b[p1] == b[i]
         nunique -= 1

      else
         p1 += 1
         b[p1] = b[i]

      end

   end  # i

   deleteat!(b,  nunique+1:n)
   #deleteat!(ii, nunique+1:n)

   return b
end # function sortunique
