export sdiag

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
