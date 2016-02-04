#=
Ipopt call
=#

n = 4
x_L = [1.0, 1.0, 1.0, 1.0]
x_U = [5.0, 5.0, 5.0, 5.0]

m = 2
g_L = [25.0, 40.0]
g_U = [2.0e19, 40.0]

function eval_f(x)
  return x[1] * x[4] * (x[1] + x[2] + x[3]) + x[3]
end

function eval_g(x, g)
  g[1] = x[1]   * x[2]   * x[3]   * x[4]
  g[2] = x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2
end

function eval_grad_f(x, grad_f)
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

prob = createProblem(n, x_L, x_U, m, g_L, g_U, 8, 10,
                     eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h)


function createProblem(
  n::Int,                     # Number of variables
  x_L::Vector{Float64},       # Variable lower bounds
  x_U::Vector{Float64},       # Variable upper bounds
  m::Int,                     # Number of constraints
  g_L::Vector{Float64},       # Constraint lower bounds
  g_U::Vector{Float64},       # Constraint upper bounds
  nele_jac::Int,              # Number of non-zeros in Jacobian
  nele_hess::Int,             # Number of non-zeros in Hessian
  eval_f,                     # Callback: objective function
  eval_g,                     # Callback: constraint evaluation
  eval_grad_f,                # Callback: objective function gradient
  eval_jac_g,                 # Callback: Jacobian evaluation
  eval_h = nothing)           # Callback: Hessian evaluation

# Set starting solution
prob.x = [1.0, 5.0, 5.0, 1.0]
addOption(prob, "hessian_approximation", "limited-memory")

# Solve
status = solveProblem(prob)

println(Ipopt.ApplicationReturnStatus[status])
println(prob.x)
println(prob.obj_val)



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

