using Ipopt

# hs071
# min x1 * x4 * (x1 + x2 + x3) + x3
# st  x1 * x2 * x3 * x4 >= 25
#     x1^2 + x2^2 + x3^2 + x4^2 = 40
#     1 <= x1, x2, x3, x4 <= 5
# Start at (1,5,5,1)
# End at (1.000..., 4.743..., 3.821..., 1.379...)

function eval_f(x) 
  return x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3]
end

function eval_g(x, g)
  # Bad: g    = zeros(2)  # Allocates new array
  # OK:  g[:] = zeros(2)  # Modifies 'in place'
  a = x[1]
  b = x[2]
  c = x[3]
  d = x[4]
  prod = a*b*c*d
  g[1] = prod
  g[2] = x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2
end

function eval_grad_f(x, grad_f)
  # Bad: grad_f    = zeros(4)  # Allocates new array
  # OK:  grad_f[:] = zeros(4)  # Modifies 'in place'
  grad_f[1] = x[1] * x[4] + x[4] * (x[1] + x[2] + x[3])
  grad_f[2] = x[1] * x[4]
  grad_f[3] = x[1] * x[4] + 1
  grad_f[4] = x[1] * (x[1] + x[2] + x[3])
end

function eval_jac_g(x, mode, rows, cols, values)
  if mode == :Structure
    # Constraint (row) 1
    rows[1] = 1; cols[1] = 1
    rows[2] = 1; cols[2] = 2
    rows[3] = 1; cols[3] = 3
    rows[4] = 1; cols[4] = 4
    # Constraint (row) 2
    rows[5] = 2; cols[5] = 1
    rows[6] = 2; cols[6] = 2
    rows[7] = 2; cols[7] = 3
    rows[8] = 2; cols[8] = 4
  else
    # Constraint (row) 1
    values[1] = x[2]*x[3]*x[4]  # 1,1
    values[2] = x[1]*x[3]*x[4]  # 1,2
    values[3] = x[1]*x[2]*x[4]  # 1,3
    values[4] = x[1]*x[2]*x[3]  # 1,4
    # Constraint (row) 2
    values[5] = 2*x[1]  # 2,1
    values[6] = 2*x[2]  # 2,2
    values[7] = 2*x[3]  # 2,3
    values[8] = 2*x[4]  # 2,4
  end
end

function eval_h(x, mode, rows, cols, obj_factor, lambda, values)
  if mode == :Structure
    # Symmetric matrix, fill the lower left triangle only
    idx = 1
    for row = 1:4
      for col = 1:row
        rows[idx] = row
        cols[idx] = col
        idx += 1
      end
    end
  else
    # Again, only lower left triangle
    # Objective
    values[1] = obj_factor * (2*x[4])  # 1,1
    values[2] = obj_factor * (  x[4])  # 2,1
    values[3] = 0                      # 2,2
    values[4] = obj_factor * (  x[4])  # 3,1
    values[5] = 0                      # 3,2
    values[6] = 0                      # 3,3
    values[7] = obj_factor * (2*x[1] + x[2] + x[3])  # 4,1
    values[8] = obj_factor * (  x[1])  # 4,2
    values[9] = obj_factor * (  x[1])  # 4,3
    values[10] = 0                     # 4,4

    # First constraint
    values[2] += lambda[1] * (x[3] * x[4])  # 2,1
    values[4] += lambda[1] * (x[2] * x[4])  # 3,1
    values[5] += lambda[1] * (x[1] * x[4])  # 3,2
    values[7] += lambda[1] * (x[2] * x[3])  # 4,1
    values[8] += lambda[1] * (x[1] * x[3])  # 4,2
    values[9] += lambda[1] * (x[1] * x[2])  # 4,3

    # Second constraint
    values[1]  += lambda[2] * 2  # 1,1
    values[3]  += lambda[2] * 2  # 2,2
    values[6]  += lambda[2] * 2  # 3,3
    values[10] += lambda[2] * 2  # 4,4
  end
end

n = 4
x_L = [1.0, 1.0, 1.0, 1.0]
x_U = [5.0, 5.0, 5.0, 5.0]

m = 2
g_L = [25.0, 40.0]
g_U = [2.0e19, 40.0]

prob = createProblem(n, x_L, x_U, m, g_L, g_U, 8, 10,
                     eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h)

prob.x = [1.0, 5.0, 5.0, 1.0]
status = solveProblem(prob)

println(Ipopt.ApplicationReturnStatus[status])
println(prob.x)
println(prob.obj_val)

#=============================================================================================================================================#


# Testing values 
beta = beta_logit;
alpha = alpha_logit;
piInc = ones(K+1,1);
piAge = ones(K+1,1);
sigma = ones(K+1,1);
xi = ones(J,1)


# Gradient of the objective function [0_theta 0_xi 2Wg]

function eval_grad(K, J, L, W, g)

	len = 1 + K + 3*(K+1) + J # = length of (theta,xi)
	up = zeros(len,1)
	down = 2*W*g 
	grad_obj = [up; down] 

end


function eval_jac(alpha, beta, piInc, piAge, sigma, x, p, inc, age, v, xi, K, N, J, L, tau)
# Defining some useful variables
	
	alpha =


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

end

## Sparsity structure ##

A = ones(J,1) #alpha
B = zeros(J,1) #beta1
C = ones(J,len_theta - 2) # rest of theta
D = ones(J,J) # xi
E = ones(L,J) # iv
F = ones(L,L) # I_g


sparse = [A B C D zero_1 ; zero_2 E  F]
