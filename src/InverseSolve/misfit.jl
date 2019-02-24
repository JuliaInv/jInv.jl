export SSDFun, HuberFun, MVFun, MVFunTotal

"""
	mis,dmis,d2mis = SSDFun(dc,dobs,Wd)

	Input:

		dc::Array   -  simulated data
		dobs::Array -  measured data
		Wd::Array   -  diagonal weighting

	Output:

		mis::Real   -  misfit, 0.5*|dc-dobs|_Wd^2
		dmis        -  gradient
		d2mis       -  diagonal of Hessian

"""
function SSDFun(dc::Union{Array{Float64},Array{Float32}},dobs::Union{Array{Float64},Array{Float32}},Wd::Union{Array{Float64},Array{Float32}})
	res   = vec(dc)-vec(dobs) # predicted - observed data
	Wd    = vec(Wd)
	mis   = .5*real(dot(Wd.*res,Wd.*res))  # data misfit
	dmis  = Wd.*(Wd.*res)
	d2mis = Wd.*Wd
	return mis, dmis, d2mis
end # function SSDFun

"""
	For complex data misfit is computed as 0.5*|real(dc)-(dobs)|_Wd^2 +  0.5*|complex(dc)-complex(dobs)|_W^2
"""
function SSDFun(dc::Array{ComplexF64},dobs::Array{ComplexF64},Wd::Array{ComplexF64})

	wdr   = vec(real(Wd)); wdi = vec(imag(Wd))
	# wdr.*dRe + im*wdi.*dIm
	res   = vec(dc)-vec(dobs)
	resw  = wdr.*real(res) + im*wdi.*imag(res)

	mis   = .5*real(dot(resw,resw))
	dmis  = sdiag(wdr.*wdr)*real(res) + 1im*sdiag(wdi.*wdi)*imag(res)
	d2mis =  wdr.*wdr + im*wdi.*wdi
	return mis,dmis,d2mis
end

"""
	mis,dmis,d2mis = HuberFun(dc,dobs,Wd,C)

	Computes misfit via

		misfit(dc,dobs) = sqrt(abs(Wd*res).^2 + eps)

	Input:
		dc::Array   -  simulated data
		dobs::Array -  measured data
		Wd::Array   -  diagional weighting
		eps         -  conditioning parameter (default=1e-3)

	Output:
		mis::Real   -  misfit
		dmis        -  gradient
		d2mis       -  diagonal of Hessian

"""
function HuberFun(dc::Array{Float64},dobs::Array{Float64},Wd::Array{Float64},eps=1e-3)
	# compute Huber distance
	res   = vec(dc-dobs)
	G     = sqrt.( abs.(Wd.*res).^2 .+ eps)
	mis   = sum(G)
	dmis  = sdiag(Wd./G)*(Wd.*res)
	d2mis = (Wd.*Wd)./G
	return mis,dmis,d2mis
end
# Old Huber fun
# function HuberFun(dc::Array{Float64},dobs::Array{Float64},Wd::Array{Float64},eps=1e-3)
# 	# compute Huber distance
# 	res   = dc-dobs
# 	G     = sqrt.( abs.(Wd.*res).^2 .+ eps)
# 	mis   = sum(G)
# 	dmis  = sdiag(Wd./G)*(Wd.*res)
# 	d2mis = (Wd.^Wd)./G
# 	return mis,dmis,d2mis
# end
