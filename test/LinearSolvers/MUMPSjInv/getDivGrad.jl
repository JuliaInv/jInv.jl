#----------------- Get the A matrix
function getDivGrad(n1,n2,n3)

        # the Divergence
        D1 = kron(sparse(I,n3,n3),kron(sparse(I,n2,n2),ddx(n1)))
        D2 = kron(sparse(I,n3,n3),kron(ddx(n2),sparse(I,n1,n1)))
        D3 = kron(ddx(n3),sparse(I,n1*n2,n1*n2))
        # DIV from faces to cell-centers
        Div = [D1 D2 D3]

        return Div*Div';
end
#----------------- 1D finite difference on staggered grid
function ddx(n)
# generate 1D derivatives
	II, JJ, VV = SparseArrays.spdiagm_internal(0 => -ones(n), 1 => ones(n)); 
	return sparse(II, JJ, VV, n, n+1)

end

