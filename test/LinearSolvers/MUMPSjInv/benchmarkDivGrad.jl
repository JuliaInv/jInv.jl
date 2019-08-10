using MUMPSjInv
using Test
using LinearAlgebra
using SparseArrays
using Printf

include("getDivGrad.jl");
Ns = (8,16,24,32,48,64,72)
MUMPStime=zeros(length(Ns))
MUMPStimeOOC=zeros(length(Ns))
Juliatime=zeros(length(Ns))


for i=1:length(Ns)
	println(@sprintf("Solve div-grad on %d x %d x %d grid",Ns[i],Ns[i],Ns[i]));
	
	A = getDivGrad(Ns[i],Ns[i],Ns[i]);
	n = size(A,1);
    rhs = randn(n);
	
	# solve using MUMPSjInv, out-of-core
	MUMPStimeOOC[i]= @elapsed begin
		x1 = solveMUMPS(A,rhs,1,1);
	end

	# solve using MUMPSjInv in core
	MUMPStime[i]= @elapsed begin
		x2 = solveMUMPS(A,rhs,1,0);
	end

	# solve using julia
	Juliatime[i] = @elapsed begin
		x3 = A\rhs;
	end
	
	@test norm(x1-x2)/norm(x1) < 1e-9
	@test norm(x1-x3)/norm(x1) < 1e-9
end
println(@sprintf("| %4s | %8s | %8s | %8s | %8s |","n","MUMPSooc","MUMPS","Julia","speedup"))

for i=1:length(Ns)
	println(@sprintf("| %d^3 | %1.5f | %1.5f | %1.5f | %1.5f |",Ns[i],MUMPStimeOOC[i],MUMPStime[i],Juliatime[i],Juliatime[i]/MUMPStime[i]))
end