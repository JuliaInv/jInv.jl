export sdiag, spoutput

import Base.kron
import Base.SparseMatrixCSC
export kron

function sdiag(a::Vector)
# S = sdiag(s) builds sparse diagonal matrix
	n = length(a)
	i = collect(1:n+1) # colptr
	j = collect(1:n)   # rowval
	return SparseMatrixCSC(n,n,i,j,a)
end

kron(v::SparseVector,A::SparseMatrixCSC) = kron(SparseMatrixCSC(v),A)
kron(A::SparseMatrixCSC,v::SparseVector) = kron(A,SparseMatrixCSC(v))

function kron(v1::SparseVector,v2::SparseVector)  
  v = kron(SparseMatrixCSC(v1),SparseMatrixCSC(v2))
   return SparseVector(v.n,v.nzind,v.nzval)
end

#--------------------------------------------------------------------------

function spoutput( filename::String,
                   A::SparseMatrixCSC )
# Output a 3 (or 4 for complex) column sparse matrix file.

f = open(filename, "w")
n = size(A,2)

complexvalue = typeof(A.nzval[1]) == Complex128

for ir = 1:n
   j1 = A.colptr[ir]
   j2 = A.colptr[ir+1] - 1

   if complexvalue
      for ic = j1:j2
         println(f, A.rowval[ic], " ", ir, " ", real(A.nzval[ic]), " ", imag(A.nzval[ic]) )
      end # ic
   else
      for ic = j1:j2
         println(f, A.rowval[ic], " ", ir, " ", A.nzval[ic] )
      end
   end # ic

end # ir

close(f)
return
end # function spoutput
