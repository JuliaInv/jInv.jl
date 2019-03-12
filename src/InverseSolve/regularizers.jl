export diffusionReg, wdiffusionReg, wTVReg,anisoTVReg, wdiffusionRegNodal,wTVRegNodal,computeRegularizer,smallnessReg,logBarrier,logBarrierSquared

function computeRegularizer(regFun::Function,mc::Vector,mref::Vector,MInv::AbstractMesh,alpha)
	R,dR,d2R = regFun(mc,mref,MInv)
	return alpha*R, alpha*dR, alpha*d2R
end

function computeRegularizer(regFun::Array{Function},mc::Vector,mref::Array,MInv::AbstractMesh,alpha::Array{Float64})
	numReg = length(regFun)
	if size(mref,2)!=numReg;
		error("computeRegularizer: number of regularizer (=$numReg) does not match number of reference models (=$(size(mref,2))).")
	end
	if length(alpha)!=numReg;
		error("computeRegularizer: number of regularizer (=$numReg) does not match number of alphas (=$(length(alpha))).")
	end

	R = 0.0; dR = zeros(length(mc)); d2R = spzeros(length(mc),length(mc))
	for k=1:numReg
		Rt,dRt,d2Rt = regFun[k](mc,mref[:,k],MInv)
		R   += alpha[k]*Rt
		dR  += alpha[k]*dRt
		d2R += alpha[k]*d2Rt
	end
	return R, dR, d2R
end


"""
	Rc,dR,d2R = diffusionReg(m,mref,M,Iact=1.0)

	Compute diffusion regularizer
		0.5*||GRAD*(m-mref)||_V^2

	Input:
		m     - model
		mref  - reference model
		M     - Mesh
		Iact  - projector on active cells

	Output
		Rc    - value of regularizer
		dR    - gradient w.r.t. m
		d2R   - Hessian
"""
function diffusionReg(m::Vector,mref,M::AbstractMesh;Iact=1.0)
	# Rc = .5* || Grad*m ||^2
	dm   = m .- mref
	Div  = getDivergenceMatrix(M)
	Div  = Iact'*Div   # project to the active cells
	V    = getVolume(M)
	Af   = getFaceAverageMatrix(M)
	d2R  = Div * sdiag(Af'*vec(Vector(diag(V)))) * Div'
	dR   = d2R*dm
	Rc   = 0.5*dot(dm,dR)
	return Rc,dR,d2R
end

"""
	Rc,dR,d2R = smallnessReg(m,mref,M,Iact=1.0)

	Compute smallness regularizer (L2 difference to reference model)

		R(m) = 0.5*||m-mref||_V^2

	Input:
		m     - model
		mref  - reference model
		M     - Mesh

	Output
		Rc    - value of regularizer
		dR    - gradient w.r.t. m
		d2R   - Hessian
"""
function smallnessReg(m::Vector,mref,M::AbstractMesh)
	dm   = m .- mref
	dR   = dm
	Rc   = 0.5*dot(dm,dm)
	return Rc,dR,sparse(1.0I,length(m),length(m))
end

function wdiffusionReg(m::Vector, mref::Vector, M::AbstractMesh; Iact=1.0, C=[])
   # Rc = a1*\\Dx(m-mref)||^2 + .. + a3*\\Dz(m-mref)||^2 + a4*|| m -mref ||^2

   if (M.dim==3) && (length(C) == 0)
      alpha1 = 1; alpha2 = 1; alpha3 = 1; alpha4 = 1e-4;
      Wt   = sdiag([alpha1*ones(M.nf[1]);alpha2*ones(M.nf[2]);alpha3*ones(M.nf[3])])
   elseif (M.dim==2) && (length(C) == 0)
	  alpha1 = 1; alpha2 = 1;alpha4 = 1e-4;
	  Wt   = sdiag([alpha1*ones(M.nf[1]);alpha2*ones(M.nf[2])])
   elseif (M.dim==3) && (length(C) == 4)
      # C = alphax, alphay, alphaz, alphas
      alpha1 = C[1]; alpha2 = C[2]; alpha3 = C[3]; alpha4 = C[4]
      Wt   = sdiag([alpha1*ones(M.nf[1]);alpha2*ones(M.nf[2]);alpha3*ones(M.nf[3])])
   elseif length(C) == 5
      #Af  = getFaceAverageMatrix(M)
      #wt  = Af'*C[5]
      alpha1 = C[1]; alpha2 = C[2]; alpha3 = C[3]; alpha4 = C[4]
      wt   = C[5]
      Wt   = sdiag([alpha1*ones(M.nf[1]);alpha2*ones(M.nf[2]);alpha3*ones(M.nf[3])]) *
             sdiag(wt)
      #Wt = sdiag(C[5])
   elseif length(C) == sum(M.nf) + 1
      # C contains alphas followed by interface weights.
      alpha4 = C[1]
      Wt = sdiag( C[2:end] )
      #println("inteface weights used.")
   else
      error("bad regparams in wdiffusionReg")
   end

   dm   = m .- mref
   Div  = getDivergenceMatrix(M)

   V    = getVolume(M)
   #Div  = Iact'*(Div*Wt)   # project to the active cells
   Div  = Iact'*(V*Div*Wt)   # project to the active cells

   Af   = getFaceAverageMatrix(M)

   #d2R  = Div * sdiag(Af'*Vector(diag(V))) * Div' + alpha4*Iact'*V*Iact
   d2R  = Div * sdiag(1 ./(Af'*Vector(diag(V)))) * Div' + alpha4*Iact'*V*Iact

   dR   = d2R*dm
   Rc   = 0.5*dot(dm,dR)

   return Rc,dR,d2R
end # function wdiffusionReg

"""
	Rc,dR,d2R = wTVReg(m,mref,M,Iact,C=[])

	Compute weighted total variation regularizer

	Input:
		m     - model
		mref  - reference model
		M     - Mesh
		Iact  - projector on active cells
		C     - anisotropy parameters (default: [1 1 1])
		eps   - conditioning parameter for TV norm (default: 1e-3)

	Output
		Rc    - value of regularizer
		dR    - gradient w.r.t. m
		d2R   - Hessian
"""
function wTVReg(m::Vector,mref,M::AbstractMesh; Iact=1.0, C=[],eps=1e-3)
	# Rc = sqrt(a1*\\Dx(m-mref)||^2 + .. + a3*\\Dz(m-mref)||^2) + a4*|| m -mref ||^2
	if length(C) == 0
		alpha1 = 1; alpha2 = 1; alpha3 = 1;
	elseif length(C) == 1
		alpha1 = 1; alpha2 = 1; alpha3 = 1;
	else
		alpha1 = C[1]; alpha2 = C[2]; alpha3 = C[3];
	end
	dm   = m .- mref
	Div  = getDivergenceMatrix(M)
	if M.dim==3
		Wt   = sdiag([alpha1*ones(M.nf[1]);alpha2*ones(M.nf[2]);alpha3*ones(M.nf[3])])
	elseif M.dim==2
		Wt   = sdiag([alpha1*ones(M.nf[1]);alpha2*ones(M.nf[2])])
	end
	Div  = Iact'*(Div*Wt)   # project to the active cells
	V    = getVolume(M); v = Vector(diag(V))
	Af   = getFaceAverageMatrix(M)

	wTV  = sqrt.(Af*(Div'*dm).^2 .+eps);

	Rc   = dot(v,wTV);
	d2R  = Div*(sdiag(Af'*(v./wTV))*Div')
	dR   = d2R*dm;
	return Rc,dR,d2R
end

"""
	Rc,dR,d2R = anisoTVReg(m,mref,M,Iact,eps)

	Compute anisotropic total variation regularizer

	Input:
		m     - model
		mref  - reference model
		M     - Mesh
		Iact  - projector on active cells
		eps   - conditioning parameter for TV norm (default: 1e-3)

	Output
		Rc    - value of regularizer
		dR    - gradient w.r.t. m
		d2R   - Hessian
"""
function anisoTVReg(m::Vector,mref,M::AbstractMesh; Iact=1.0, eps=1e-3)
	# Rc = sqrt(a1*\\Dx(m-mref)||^2 + .. + a3*\\Dz(m-mref)||^2) + a4*|| m -mref ||^2
	dm   = m .- mref
	Div  = getVolume(M)*getDivergenceMatrix(M)
	Div  = Iact'*(Div)   # project to the active cells
	v    = ones(size(Div,2)) 
	
	wTV  = sqrt.((Div'*dm).^2 .+eps);
	Rc   = dot(v,wTV);
	d2R  = Div*(sdiag(v./wTV)*Div')
	dR   = d2R*dm;
	return Rc,dR,d2R
end

"""
	Rc,dR,d2R = wdiffusionRegNodal(m::Vector, mref::Vector, M::AbstractMesh; Iact=1.0, C=[])

	Computes weighted diffusion regularizer for nodal model

	Input:
		m     - model
		mref  - reference model
		M     - Mesh
		Iact  - projector on active cells
		C     - optional parameters

	Output
		Rc    - value of regularizer
		dR    - gradient w.r.t. m
		d2R   - Hessian
"""
function wdiffusionRegNodal(m::Vector, mref::Vector, M::AbstractMesh; Iact=1.0, C=[])
	dm = m.-mref;
	if isempty(C)
		C = [1,1,1,1,1e-5];
	end
	if M.dim==3
		Wt   = [C[1]*ones(M.ne[1]);C[2]*ones(M.ne[2]);C[3]*ones(M.ne[3])];
	else
		Wt   = [C[1]*ones(M.ne[1]);C[3]*ones(M.ne[2])];
	end

	V    = getVolume(M);
	Af   = getEdgeAverageMatrix(M)
	v    = (Af'*Vector(diag(V)))
	Wt   = sdiag(Wt.*(v.*Wt));

	Av = getNodalAverageMatrix(M);

	Av = Av*Iact;
	mass = Av'*V*Av;

	Grad   = getNodalGradientMatrix(M)*Iact;

	d2R = Grad'*Wt*Grad;
	d2R += C[4]*mass;
	dR  = d2R*dm;
	Rc  = 0.5*dot(dm,dR);

	if isnan(Rc)
		dump(m)
	end
	clear!(M);
   return Rc,dR,d2R
end


function wTVRegNodal(m::Vector, mref::Vector, M::AbstractMesh; Iact=1.0, C=[])
	dm = m.-mref;
	if isempty(C)
		C = [1,1,1,1,1e-5];
	end
	eps = 1e-3;
	if M.dim==3
		Wt   = sdiag([C[1]*ones(M.ne[1]);C[2]*ones(M.ne[2]);C[3]*ones(M.ne[3])]);
	else
		Wt   = sdiag([C[1]*ones(M.ne[1]);C[3]*ones(M.ne[2])]);
	end


	V    = getVolume(M);
	Av = getNodalAverageMatrix(M);
	mass = Iact'*Av'*V*Av*Iact;
	d2R  = C[4]*mass;
	Rc   = 0.5*dot(dm,d2R*dm);

	v      = Vector(diag(V))
	Af     = getEdgeAverageMatrix(M)
	Grad   = Wt*getNodalGradientMatrix(M)*Iact;
	wTV    = sqrt.(Af*(((Grad*dm).^2).+eps));
	Rc     += dot(v,wTV);
	d2R    += Grad'*sdiag(Af'*(v./wTV))*Grad;
	dR     = d2R*dm
	if isnan(Rc)
		dump(m)
	end
	clear!(M);
   return Rc,dR,d2R
end


"""
	Rc,dR,d2R = logBarrier(m::Vector, z::Vector, M::AbstractMesh,low::Vector,high::Vector, epsilon)

	Computes logBarrier regularizer

	R = -log(1 - ((m-high)/epsilon).^2) if high-epsilon <  m < high
					0					if low+epsilon  <= m <= high-epsilon
		-log(1 - ((m-low)/epsilon).^2)	if low          <  m < low+epsilon

	Input:
		m     	- model
		z     	- not being used. Here for compatibility.
		M     	- Mesh. not being used. Here for compatibility.
		low   	- low bound for each coordinate.
		high  	- high bound for each coordinate.
		epsilon - layer width of the barier.

	Output
		g    - value of regularizer
		dg    - gradient w.r.t. m
		d2g   - Hessian (diagonal matrix). Second derivative is not continous.
"""
function logBarrier(m::Vector,z::Vector,M::AbstractMesh, low::Vector, high::Vector ,epsilon = min(0.1*abs(low),0.1*abs(high)))
	
	if length(findall(m.>=high)) + length(findall(m.<=low))>0
		return Inf, zeros(length(m)), spzeros(length(m),length(m));
	end
	
	z = copy(m);
	indProj = zeros(length(m));

	low = low + epsilon;
	high = high - epsilon;
	t = z.<low;
	indProj[t]  .= 1.0;
	z[t]  = low[t];
	t = z.>high;
	indProj[t] .= 1.0;
	z[t] = high[t];

	e 	= (epsilon.+1e-10);
	dm  = (m-z)./e;
	f   = 1.0 .- dm.^2; # > 0
	if minimum(f) < 0.0
		error("negative f");
	end
	df  = (-2.0).*indProj.*(dm./e);
	d2f = (-2.0).*(indProj./(e.^2));

	g   = -sum(log.(f));
	dg  = -df./f;
	d2g = dg.^2 - d2f./f  ;
	return g,dg,sdiag(d2g);
end


"""
	Rc,dR,d2R = logBarrierSquared(m::Vector, z::Vector, M::AbstractMesh,low::Vector,high::Vector, epsilon)

	Computes logBarrier regularizer

	R = (log(1 - ((m-high)/epsilon).^2))^2 if high-epsilon <  m < high
					0					   if low+epsilon  <= m <= high-epsilon
		(log(1 - ((m-low)/epsilon).^2))^2  if low          <  m < low+epsilon

	Input:
		m     	- model
		z     	- not being used. Here for compatibility.
		M     	- Mesh. not being used. Here for compatibility.
		low   	- low bound for each coordinate.
		high  	- high bound for each coordinate.
		epsilon - layer width of the barier.

	Output
		g    - value of regularizer
		dg    - gradient w.r.t. m
		d2g   - Gauss Newton Hessian approximation (diagonal matrix). Second derivative approx is continous.
"""
function logBarrierSquared(m::Vector,z::Vector,M::AbstractMesh, low::Vector, high::Vector ,epsilon = min(0.1*abs(low),0.1*abs(high)))
	if length(findall(m.>=high)) + length(findall(m.<=low))>0
		return Inf, zeros(length(m)), spzeros(length(m),length(m));
	end
	z = copy(m);
	indProj = zeros(length(m));
	low = low + epsilon;
	high = high - epsilon;
	t = z.<low;
	indProj[t]  .= 1.0;
	z[t]  = low[t];
	t = z.>high;
	indProj[t] .= 1.0;
	z[t] = high[t];

	e 	= (epsilon.+1e-10);
	dm  = (m-z)./e;
	f   = 1.0 .- dm.^2; # > 0
	if minimum(f) < 0.0
		error("negative f");
	end
	df  = (-2.0).*indProj.*(dm./e);

	r   = log.(f);
	dr  = df./f;
	g   = 0.5*dot(r,r);
	dg  = dr.*r;
	d2g = dr.^2;
	return g,dg,sdiag(d2g);
end

export TikhonovReg
function TikhonovReg(m::Vector,mref::Vector,M::AbstractMesh,sigma::SparseMatrixCSC)
d = (m-mref);
df = sigma*d;
f = 0.5*dot(d,df);
d2f = sigma;
return f,df,d2f;
end
