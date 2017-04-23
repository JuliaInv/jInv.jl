
export sortpermFast


function sortpermFast(A::Vector)
   const n = length(A)

   ii = collect(1:n)
   B = copy(A)
   quicksort!(B,ii, 1,n)

   return ii, B  #  B = A[ii]
end # function sortpermFast

#----------------------------------------------------

function sortpermFast(A::Vector, D::Vector)
   # Sort A and permute D according to A.
   # For duplicate values in A, keep only values corresponding to
   # the SMALLEST D.
   const n = length(A)
   if length(D) != n
      error("Lengths of A and D must be the same.")
   end

   ii = collect(1:n)
   quicksort!(A, ii, 1,n)

   D = D[ii]
   
   if allunique(A)
      return A, D
   end
   
   idxkeep = trues(n)
   
   for j = 2 : n
      for k = j-1 : -1 : 1
         if A[j] != A[k]
            break
         end
         
         if D[j] < D[k]
            idxkeep[k] = false
         else
            idxkeep[j] = false
         end
         
      end  # k
   end  # j

   A = A[idxkeep]
   D = D[idxkeep]

   return A, D
end # function sortpermFast

#----------------------------------------------------

function quicksort!(A, order, i=1,j=length(A))
# modified from:
# http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Julia	

	 @inbounds begin
    if j > i
    	
    	  if  j - i <= 10 
    	  	  # Insertion sort for small groups is faster than Quicksort
    	     InsertionSort!(A,order, i,j)
    	     return A
    	  end
    	
        #pivot = A[rand(i:j)] # random element of A
        pivot = A[ div(i+j,2) ] 
        left, right = i, j
        while left <= right
            while A[left] < pivot
                left += 1
            end
            while A[right] > pivot
                right -= 1
            end
            if left <= right
                A[left], A[right] = A[right], A[left]
                order[left], order[right] = order[right], order[left]

                left += 1
                right -= 1
            end
        end  # left <= right
        
        quicksort!(A,order, i,   right)
        quicksort!(A,order, left,j)
    end  # j > i
    end
    
    return A
end # function quicksort!

#----------------------------------------------------

function InsertionSort!(A, order, ii=1, jj=length(A))
	
	 @inbounds begin
    for i = ii+1 : jj
        j = i - 1
        temp  = A[i]
        itemp = order[i]
        
        while true
            if j == ii-1
            	 break
            end
            if A[j] <= temp
            	 break
            end
            A[j+1] = A[j]
            order[j+1] = order[j]
            j -= 1
        end
        
        A[j+1] = temp
        order[j+1] = itemp
    end  # i
    end

return
end # function InsertionSort!
