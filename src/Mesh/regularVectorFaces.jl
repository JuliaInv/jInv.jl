export getUjProjMatrix
export getVectorMassMatrix, getDifferentialOperators,GetLinearElasticityOperator,GetLinearElasticityOperatorFullStrain

function ddxCN(n::Int64,h) ## From nodes to cells
	D = (1/h)*ddx(n);
	return D
end

function ddxNC(n::Int64,h) ## From cells to nodes
	#A = 1/h*spdiagm((-ones(n),ones(n)),[-1,0],n+1,n)
	I, J, V = SparseArrays.spdiagm_internal(-1 => fill(-1.0,n), 0 => fill(1.0,n)); 
	D = sparse(I, J, (1/h)*V,n+1,n);
	## This is for boundary conditions: need to figure this out
	D[1,1] = 1.0/h;
	D[end,end] = -1.0/h;
	return D
end

 	 
function getDifferentialOperators(M::RegularMesh,optype = 1)
	n = M.n;
	h = M.h;
	if length(n)==3
		# Face sizes
		nf1 = prod(n + [1; 0; 0])
		nf2 = prod(n + [0; 1; 0])
		nf3 = prod(n + [0; 0; 1])
		nf  = [nf1; nf2; nf3];
		
		# Notation Dij = derivative of component j in direction i
		tmp = ddxCN(n[1],h[1])
		D11 = kron(speye(n[3]),kron(speye(n[2]),tmp))

		tmp = ddxNC(n[1],h[1])
		D12 = kron(speye(n[3]),kron(speye(n[2]+1),tmp))

		tmp = ddxNC(n[1],h[1])
		D13 = kron(speye(n[3]+1),kron(speye(n[2]),tmp))

		tmp = ddxNC(n[2],h[2])
		D21 = kron(speye(n[3]),kron(tmp,speye(n[1]+1)))

		tmp = ddxCN(n[2],h[2])
		D22 = kron(speye(n[3]),kron(tmp,speye(n[1])))

		tmp = ddxNC(n[2],h[2])
		D23 = kron(speye(n[3]+1),kron(tmp,speye(n[1])))

		tmp = ddxNC(n[3],h[3])
		D31 = kron(tmp,kron(speye(n[2]),speye(n[1]+1)))

		tmp = ddxNC(n[3],h[3])
		D32 = kron(tmp,kron(speye(n[2]+1),speye(n[1])))

		tmp = ddxCN(n[3],h[3])
		D33 = kron(tmp,kron(speye(n[2]),speye(n[1])))

		# Tensor sizes
		t = [size(D11,1); size(D12,1); size(D13,1);
			size(D21,1); size(D22,1); size(D23,1);
			size(D31,1); size(D32,1); size(D33,1);]
		
		
		Div = [D11 D22 D33]
		if optype==1
			vectorDiv = [Div; spzeros(t[2]+t[3]+t[4],sum(nf)); Div; spzeros(t[6]+t[7]+t[8],sum(nf)); Div];
			vectorGrad  = 0.5*[blockdiag(D11,D12,D13);blockdiag(D21,D22,D23);blockdiag(D31,D32,D33)]
			vectorGradT = 0.5*blockdiag([D11;D21;D31],[D12;D22;D32],[D13;D23;D33])
			D11 = 0; D12 = 0; D13 = 0; D22 = 0; D21 = 0; D23 = 0; D31 = 0; D32 = 0; D33 = 0;
			vectorGrad  = (vectorGrad  + vectorGradT);
		else
			vectorGrad = blockdiag([D11;D21;D31],[D12;D22;D32],[D13;D23;D33]);
			vectorDiv = 0;
			vectorGradT = 0;
		end
	else		
		nf1 = prod(n + [1; 0]);
		nf2 = prod(n + [0; 1]);
		nf  = [nf1; nf2];
		
		# Notation Dij = derivative of component j in direction i
		tmp = ddxCN(n[1],h[1]);
		D11 = kron(speye(n[2]),tmp);

		tmp = ddxNC(n[1],h[1]);
		D12 = kron(speye(n[2]+1),tmp);

		tmp = ddxNC(n[2],h[2])
		D21 = kron(tmp,speye(n[1]+1))

		tmp = ddxCN(n[2],h[2])
		D22 = kron(tmp,speye(n[1]))

		# Tensor sizes
		t = [size(D11,1); size(D12,1); size(D21,1); size(D22,1);]
		  
		Div = [D11 D22]
		
		
		if optype==1
			vectorDiv = [Div; spzeros(t[2]+t[3],sum(nf)); Div; ];
			vectorGrad  = [blockdiag(D11,D12);blockdiag(D21,D22)]
			vectorGradT = blockdiag([D11;D21],[D12;D22])
			vectorGrad  = 0.5*(vectorGrad  + vectorGradT)
		else
			# vectorGrad = [[D11;D21],[D12;D22]];
			vectorGrad = blockdiag([D11;D21],[D12;D22]); ## diagGrad
			vectorDiv = 0;
			vectorGradT = 0;
		end
	end 
	return vectorGrad, vectorDiv, Div,nf,vectorGradT
end






function getTensorMassMatrix(M::RegularMesh, mu::Vector; saveAvMat::Bool = true)
 
	if M.dim==3
		n = M.n;
		# (Ae1, Ae2, Ae3) = getEdgeAveragingMatrices(n)
		# mu    = vec(mu);
		# Ae3mu = Ae3*mu;
		# Ae1mu = Ae1*mu;
		# Ae2mu = Ae2*mu;
		
		Ae = getEdgeAverageMatrix(M; saveMat = saveAvMat, avN2C = avN2C_Nearest);
		AeMu = Ae'*mu;
		Ae1mu = AeMu[1:M.nf[1]];
		Ae2mu = AeMu[(M.nf[1]+1):(M.nf[1]+M.nf[2])];
		Ae3mu = AeMu[(M.nf[1]+M.nf[2]+1):end];		
		m = [mu;Ae3mu;Ae2mu;Ae3mu;mu;Ae1mu;Ae2mu;Ae1mu;mu];
		Mass = spdiagm(m);
		len = [0;length(mu)+length(Ae3mu)+length(Ae2mu);length(Ae3mu)+length(mu)+length(Ae1mu);length(Ae2mu)+length(Ae1mu)+length(mu)];
		len[3] = len[3]+len[2];
		len[4] = len[4] + len[3];
	else
		# (Ae1, Ae2) = getEdgeAveragingMatrices(n)
		# mu    = vec(mu);
		# Ae1mu = Ae1*mu;
		# Ae2mu = Ae2*mu;
		# m = [mu;Ae1mu;Ae2mu;mu];
		# Mass = spdiagm(m);
		# len = [0;length(mu)+length(Ae1mu);length(Ae2mu)+length(mu)];
		# len[3] = len[3]+len[2];
		
		
		# In 2D the strain tensor (when u is on faces) is defined on the nodes.
		An = getNodalAverageMatrix(M; saveMat = saveAvMat, avN2C = avN2C_Nearest)
		AnMu = An'*mu;
		Mass = spdiagm([mu;AnMu;AnMu;mu]);
		len_nodes = (M.n[1]+1)*(M.n[2]+1);
		len = [0;length(mu)+len_nodes;len_nodes+length(mu)];
		len[3] = len[3]+len[2];
	end
	return Mass,len
end


function getUjProjMatrix(n,j)
	if length(n)==3
		nf1 = prod(n + [1; 0; 0])
		nf2 = prod(n + [0; 1; 0])
		nf3 = prod(n + [0; 0; 1])
		nf  = [0;nf1;nf2+nf1;nf3+nf1+nf2];
		Pj = speye(nf1+nf2+nf3);
		Pj = Pj[:,nf[j]+1:nf[j+1]];
	else
		nf1 = prod(n + [1; 0])
		nf2 = prod(n + [0; 1])
		nf  = [0;nf1;nf2+nf1];
		Pj = speye(nf1+nf2);
		Pj = Pj[:,nf[j]+1:nf[j+1]];
	end
	return Pj;
end


function GetLinearElasticityOperator(M::RegularMesh, mu::Vector,lambda::Vector)
	vecGRAD,~,Div,nf, = getDifferentialOperators(M,2);
	Mu     = getTensorMassMatrix(M,mu,saveAvMat=false)[1];
	LambdaMuCells = spdiagm(lambda[:] + mu[:]);
	Rho	   = getFaceMassMatrix(M,1e-5*mu, saveMat = false, avN2C = avN2C_Nearest);
	G = vecGRAD'*Mu*vecGRAD
	
	H = Div'*LambdaMuCells*Div + vecGRAD'*Mu*vecGRAD + Rho;
	return H;						
end

function GetLinearElasticityOperatorFullStrain(M::RegularMesh, mu::Vector,lambda::Vector)
	vecGRAD,vecDiv,Div,nf, = getDifferentialOperators(M,1);
	Mu     = getTensorMassMatrix(M,mu,saveAvMat=false)[1];
	Lambda = getTensorMassMatrix(M,lambda,saveAvMat=false)[1];
	Rho	   = getFaceMassMatrix(M,1e-5*mu, saveMat = false, avN2C = avN2C_Nearest);
	H = vecGRAD'*(Lambda*vecDiv + 2.0*Mu*vecGRAD) + Rho;
	return H;						
end
