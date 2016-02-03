
# Testing values 
beta = ones(K,1);
alpha = 1.0;
piInc = ones(K,1);
piAge = ones(K,1);
sigma = ones(K,1):
xi = ones(J,1)



# Gradient of the objective function [0_theta 0_xi 2Wg]

function eval_grad(K, J, L, W, g)

	len = 1 + 4*K + J # = length of (theta,xi)
	up = zeros(len,1)
	down = 2*W*g 
	grad_obj = [up; down] 

end

# Defining some useful variables

denom = zeros(N,1);
tau = zeros(J,N);

for n = 1:N

	denom[n] = sum( exp(beta[1]
            -(alpha + piInc[1]*inc[n] + piAge[1]*age[n] + sigma[1]*v[n,1])*p
            + x[:,2:K]*(beta[2:K] + piInc[2:K]*inc[n] + piAge[2:K]*age[n] + Diagonal(sigma[2:K])*v'[2:K,n] )
            + xi )
            ) 
																				
	tau[:,n] = exp(beta[1]																# defining the tau as in the MPEC paper appendix
            -(alpha + piInc[1]*inc[n] + piAge[1]*age[n] + sigma[1]*v[n,1])*p
            + x[:,2:K]*(beta[2:K] + piInc[2:K]*inc[n] + piAge[2:K]*age[n] + Diagonal(sigma[2:K])*v'[2:K,n] )
            + xi )/denom[n]
end

	# Jacobian of the constraints [ds/dtheta  ds/dxi  0 \\ 0 -Z' eye(g)] , theta = (alpha, beta, piInc, piAge, sigma)

function eval_jac(alpha, beta, piInc, piAge, sigma, x, p, inc, age, v, xi, K, N, J, L, tau)

	# allocate matrices

	d_alpha = zeros(J,N);
	d_beta = zeros(J,K,N);
	d_piInc = zeros(J,K,N);
	d_piAge = zeros(J,K,N);
	d_sigma = zeros(J,K,N);
	d_xi = zeros(J,J,N);

	# define derivatives point by point

	for n = 1:N
		for j = 1:J
			d_alpha[j,n] = p[j]*tau[j,n] - sum(Diagonal(p)*tau[:,n])*tau[j,n]
			d_piInc[j,1,N] = (p[j]*tau[j,n] - sum(Diagonal(p)*tau[:,n])*tau[j,n])*inc[n]
			d_piAge[j,1,N] = (p[j]*tau[j,n] - sum(Diagonal(p)*tau[:,n])*tau[j,n])*age[n]
			d_sigma[j,1,N] = (p[j]*tau[j,n] - sum(Diagonal(p)*tau[:,n])*tau[j,n])*v'[1,n]
			for k=2:K
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

	D_theta = [D_alpha D_beta D_piInc D_piAge D_sigma]

	D_xi = (1/N)*sum(d_xi, 3)

	# Now that we have derivatives, the rest is easy.

	len_theta = 1 + 4*K 

	zero_1 = zeros(J,L)
	zero_2 = zeros(L,len_theta)
	I_g = eye(L)

	jac = [D_theta D_xi zero_1 ; zero_2 -iv' I_g]
end

## Sparsity structure ##

A = ones(J,1) #alpha
B = zeros(J,1) #beta1
C = ones(J,4*K-1) # rest of theta
D = ones(J,J) # xi
E = ones(L,J) # iv
F = ones(L,L) # I_g

sparse = [ A B C D zero_1 ; zero_2 E  F]