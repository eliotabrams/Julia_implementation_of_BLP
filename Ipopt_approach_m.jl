function eval_jac_g(param, mode, rows, cols, values)

# Redefine variables
	alpha = param[1];
	beta = param[2:5];
	piInc = param[6:10];
	piAge = param[11:15];
	sigma = param[16:20];
	g = param[21:30];
	xi = param[31:558];

# Defining some useful variables
	d_alpha = zeros(J,N);
	d_beta = zeros(J,K,N);
	d_piInc = zeros(J,K+1,N);
	d_piAge = zeros(J,K+1,N);
	d_sigma = zeros(J,K+1,N);
	d_xi = zeros(J,J,N);

	denom = zeros(N,1);
	tau = zeros(J,N);

	for n = 1:N

		denom[n] = sum( exp( 
	            -(alpha + piInc[K+1]*inc[n] + piAge[K+1]*age[n] + sigma[K+1]*v[n,K+1])*p
	            + x[:,1:K]*(beta[1:K] + piInc[1:K]*inc[n] + piAge[1:K]*age[n] + Diagonal(sigma[1:K])*v'[1:K,n] )
	            + xi )
	            ) 
																					
		tau[:,n] = exp(																# defining the tau as in the MPEC paper appendix
	            -(alpha + piInc[K+1]*inc[n] + piAge[K+1]*age[n] + sigma[K+1]*v[n,K+1])*p
	            + x[:,1:K]*(beta[1:K] + piInc[1:K]*inc[n] + piAge[1:K]*age[n] + Diagonal(sigma[1:K])*v'[1:K,n] )
	            + xi )/denom[n]
	end

	# Jacobian of the constraints [ds/dtheta  ds/dxi  0 \\ 0 -Z' eye(g)] , theta = (alpha, beta, piInc, piAge, sigma

	# allocate matrices

	# define derivatives point by point

	for n = 1:N
		for j = 1:J
			d_alpha[j,n] = p[j]*tau[j,n] - sum(Diagonal(p)*tau[:,n])*tau[j,n]
			d_piInc[j,K+1,N] = (p[j]*tau[j,n] - sum(Diagonal(p)*tau[:,n])*tau[j,n])*inc[n]
			d_piAge[j,K+1,N] = (p[j]*tau[j,n] - sum(Diagonal(p)*tau[:,n])*tau[j,n])*age[n]
			d_sigma[j,K+1,N] = (p[j]*tau[j,n] - sum(Diagonal(p)*tau[:,n])*tau[j,n])*v'[K+1,n]
			for k=1:K
				d_beta[j,k,n] = ( x[j,k]*tau[j,n] - sum(Diagonal(x[:,k])*tau[:,n])*tau[j,n] ) # this is awesome!
				d_piInc[j,k,n] = ( x[j,k]*tau[j,n] - sum(Diagonal(x[:,k])*tau[:,n])*tau[j,n] )*inc[n]
				d_piAge[j,k,n] = ( x[j,k]*tau[j,n] - sum(Diagonal(x[:,k])*tau[:,n])*tau[j,n] )*age[n]
				d_sigma[j,k,n] = ( x[j,k]*tau[j,n] - sum(Diagonal(x[:,k])*tau[:,n])*tau[j,n] )*v'[k,n]
			end
			for jj=1:J
				if j == jj
					d_xi[j,jj,n] = tau[j,n]*(1 - tau[j,n])
				else 
					d_xi[j,jj,N] = - tau[j,n]*tau[jj,n]
				end 
			end
		end
	end

	# define derivative matrices

	D_alpha = (1/N)*sum(d_alpha, 2)
	D_beta = (1/N)*sum(d_beta, 3) #note d_beta1 = 0
	D_piInc = (1/N)*sum(d_piInc, 3)
	D_piAge = (1/N)*sum(d_piAge, 3)
	D_sigma = (1/N)*sum(d_sigma, 3)

	D_theta = [D_beta D_alpha D_piInc D_piAge D_sigma]

	D_xi = (1/N)*sum(d_xi, 3)

	# Now that we have derivatives, the rest is easy.

	len_theta = 1 + K + 3*(K+1) 

	zero_1 = zeros(J,L)
	zero_2 = zeros(L,len_theta)
	I_g = eye(L)

	jac = [D_theta D_xi zero_1 ; zero_2 -iv' I_g]

	jac = convert(Array{Float64,2}, jac[1:size(jac,1), 1:size(jac,2)])
	(I, J, V) = findnz(jac);


if mode == :Structure
    rows = I; cols = J;
  else
    values = V;
  end
end

