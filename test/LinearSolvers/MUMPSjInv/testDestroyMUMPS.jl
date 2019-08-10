# Factorize two different matrices (one real, one complex) and solve for multiple rhs and free memory
using Distributed
using LinearAlgebra
using SparseArrays
using Test

include("getDivGrad.jl");

Ar = getDivGrad(30,31,33);
n  = size(Ar,1);
Ac  = Ar + im*spdiagm(0=>rand(n));

nrhs = 10;
rhsr = randn(n,nrhs) + im*randn(n,nrhs);
rhsc = randn(n,nrhs);


# println("Factorize real matrix on main instance and destroy")
# tic();
# fact = 0;
# for i=1:100;
# 	fact = factorMUMPS(Ar,1)
# 	destroyMUMPS(fact)
# end
# toc();

println("Factorize real matrix on remote worker and destroy")
@time begin
	fact = 0;
	for i=1:10;
		fact = remotecall_fetch(factorMUMPS,workers()[1], Ar,1)
		remotecall_fetch(destroyMUMPS,workers()[1],fact)
	end
end
println("DONE!")


