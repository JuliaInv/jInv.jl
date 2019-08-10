using MUMPSjInv
using Test
using LinearAlgebra
using SparseArrays
using Printf


include("getDivGrad.jl");
Ns = (8,16,24,32,48,64)
MUMPStime=zeros(length(Ns))
MUMPStimeOOC=zeros(length(Ns))
Juliatime=zeros(length(Ns))


for i=1:length(Ns)
	println(@sprintf("Solve div-grad on %d x %d x %d grid",Ns[i],Ns[i],Ns[i]));
	
	A = getDivGrad(Ns[i],Ns[i],Ns[i]);
	n = size(A,1);
    rhs = randn(n);
	
	# solve using MUMPSjInv
	MUMPStimeOOC[i] = @elapsed begin
		x = solveMUMPS(A,rhs,1,1);
	end	
	# solve using MUMPSjInv in core
	MUMPStime[i] = @elapsed begin
		x = solveMUMPS(A,rhs,1,0);
	end

	# solve using MUMPSjInv
	Juliatime[i] = @elapsed begin
		x = A\rhs;
	end
end

for i=1:length(Ns)
	println(@sprintf("| %d^3 | %1.3f | %1.3f | %1.3f | %1.1f |",Ns[i],MUMPStimeOOC[i],MUMPStime[i],Juliatime[i],Juliatime[i]/MUMPStime[i]))
end