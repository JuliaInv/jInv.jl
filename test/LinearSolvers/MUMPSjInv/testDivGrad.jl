# Test for solving div-grad system both complex and real and single and multiple rhs
using MUMPSjInv
using Test
using LinearAlgebra
using SparseArrays

include("getDivGrad.jl");

A = getDivGrad(32,32,16);
n = size(A,1);

# REAL: test for single rhs
println("Test for real SPD matrix: one rhs");
rhs = randn(n);

x = solveMUMPS(A,rhs,1);
err  =  norm(A*x-rhs) / norm(rhs);
@test err < 1e-14

# REAL: test for multiple rhs
println("Test for real SPD matrix: multiple rhs");
nrhs = 10;
rhs = randn(n,nrhs);

x = solveMUMPS(A,rhs,1);

err = zeros(nrhs)
for i=1:nrhs
        err[i] =  norm(A*x[:,i]-rhs[:,i]) / norm(rhs[:,i]);
end
@test eltype(x) == Float64
@test maximum(err) < 1e-14

# COMPLEX: test for multiple rhs
println("Test for complex symmetric matrix: one rhs");
r = rand(n);
A = A + im*spdiagm(0=>r);

rhs = randn(n) + im*randn(n);

x = solveMUMPS(A,rhs,1);

err=  norm(A*x-rhs) / norm(rhs);
@test eltype(x) == ComplexF64
@test err < 1e-14

# COMPLEX : test for multiple rhs
println("Test for complex symmetric matrix: multiple rhs");
nrhs = 10;
rhs = randn(n,nrhs) + im*randn(n,nrhs);

x = solveMUMPS(A,rhs,2);
err = zeros(nrhs)
for i=1:nrhs
        err[i] =  norm(A*x[:,i]-rhs[:,i]) / norm(rhs[:,i]);
end
@test eltype(x) == ComplexF64
@test maximum(err) < 1e-14

println("DONE!")

