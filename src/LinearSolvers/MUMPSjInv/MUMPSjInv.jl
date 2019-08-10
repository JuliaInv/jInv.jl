__precompile__()

module MUMPSjInv

using LinearAlgebra
using SparseArrays
using Distributed
using Printf

mutable struct MUMPSfactorization{T}
	ptr::Int64     # pointer to factorization
	worker::Int64  # id of worker that holds factorization
	n::Int64       # matrix size
	a11::T         # first element (HACK for generating the right parametric type)
	time::Float64  # factorization time
end
	const MUMPSlibPath  = abspath(joinpath(splitdir(Base.source_path())[1],"..","lib","MUMPS"))

	include("MUMPSfuncs.jl")

	export solveMUMPS,solveMUMPS!, factorMUMPS, applyMUMPS,applyMUMPS!,destroyMUMPS, MUMPSfactorization

end
