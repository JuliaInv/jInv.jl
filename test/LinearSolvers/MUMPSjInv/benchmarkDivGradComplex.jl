using MUMPSjInv
if VERSION >= v"0.5.0-dev+7720"
    using Test
else
    using BaseTestNext
    const Test = BaseTestNext
end
include("getDivGrad.jl");
Ns = (8,16,24,32,48,64)
MUMPStime=zeros(length(Ns))
Juliatime=zeros(length(Ns))
nrhs = 1

for i=1:length(Ns)
	println(@sprintf("Solve div-grad on %d x %d x %d grid",Ns[i],Ns[i],Ns[i]));
	A = getDivGrad(Ns[i],Ns[i],Ns[i]);
	n = size(A,1);
    r = ones(n);
	A = A + im*spdiagm(r,0);
	rhs = randn(n,nrhs)+ im*randn(n,nrhs)

	# solve using MUMPSjInv
	@elapsed MUMPStime[i] = begin
	    x = solveMUMPS(A,rhs,1,0);
    end

	# solve using MUMPSjInv
	Juliatime[i]= @elapsed
        x = A\rhs;
    end
end

println(@sprintf("| %4s | %8s | %8s | %8s |","n","MUMPS","Julia","speedup"))

for i=1:length(Ns)
	println(@sprintf("| %d^3  | %1.5f | %1.5f | %1.5f |",Ns[i],MUMPStime[i],Juliatime[i],Juliatime[i]/MUMPStime[i]))
end
