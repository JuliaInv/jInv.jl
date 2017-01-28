
export sortpermFast


function sortpermFast(A)
   const n = length(A)

   ii = collect(1:n)
   B = copy(A)
   quicksort!(B,ii, 1,n)

   return ii, B
end # function mysortperm

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
