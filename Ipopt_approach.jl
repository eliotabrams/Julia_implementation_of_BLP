#=
Ipopt call
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
=#

#####################
##     Setup       ##
#####################

using Ipopt
using JuMP
using DataFrames
#cd("/Users/eliotabrams/Desktop/Advanced\ Industrial\ Organization\ 2/Julia_implementation_of_BLP")


#####################
##      Data       ##
#####################

# Load data
product = DataFrames.readtable("dataset_cleaned.csv", separator = ',', header = true);
population = DataFrames.readtable("population_data.csv", separator = ',', header = true);

# Define variables
x = product[:,3:6];
p = product[:,7];
z = product[:,8:13];
s0 = product[:,14];
s = product[:,2];
iv = [x z];
inc = population[:,1];
age = population[:,2];
v = population[:,3:7];

# Store dimensions
K = size(x,2);
L = K+size(z,2);
J = size(x,1);
N = size(v,1);
M = size(v,2);


##########################
##  Simple Logit Model  ##
##########################

# Setup the simple logit model
logit = Model(solver = IpoptSolver(tol = 1e-8, max_iter = 1000, output_file = "logit2.txt"));

# Define variables
@defVar(logit, g[1:L]);
@defVar(logit, xi[1:J]);
@defVar(logit, alpha);
@defVar(logit, beta[1:K]);

# We minimize the gmm objective with the identity as the weighting matrix
# subject to the constraints g = sum_j xi_j iv_j and market share equations
@setObjective(logit, Min, sum{g[l]^2,l=1:L});
@addConstraint(
    logit, 
    constr[l=1:L], 
    g[l]==sum{xi[j]*iv[j,l], j=1:J}
);
@addNLConstraint(
    logit, 
    constr[j=1:J], 
    xi[j]==log(s[j])-log(s0[j])+alpha*p[j]-sum{beta[k]*x[j,k],k=1:K}
);

# Solve the model
status = solve(logit);

# Print the results
print(status)
print("alpha = ", getValue(alpha))
print("beta = ", getValue(beta[1:K]))

# Save results to use in the setup of BLP Model
g_logit=getValue(g);
xi_logit=getValue(xi);
alpha_logit=getValue(alpha);
beta_logit=getValue(beta);


##################################
##   Define Ipopt call for BLP  ##
##################################

# There are 558 variables that we do not know much info on
n = 558;
x_L = -1000*ones(n);
x_U = 1000*ones(n);

# We have 538 equality constraints
m = 538;
g_L = zeros(m);
g_U = zeros(m);

# Number of non-zeros
nele_jac = 294196;
nele_hess = 0; 

iv = convert(Array, iv);
W = inv((1/J)*iv'*Diagonal(diag(xi_logit*xi_logit'))*iv);

function eval_f(param)
  return param[21:30]'W*param[21:30]
end

function eval_grad_f(param, grad_f)
  grad_f = 2*W*param[21:30]
end

#=
g = [share, g]
g = [1:528, 529:538]
param = [alpha, beta,  piInc, piAge, sigma, g,     xi     ]
param = [1,     2:5,   6:10,  11:15, 16:20, 21:30, 31:558
=#
function eval_g(param, g)
  alpha = param[1];
  beta = param[2:5];
  piInc = param[6:10];
  piAge = param[11:15];
  sigma = param[16:20];
  denom = zeros(N,1);
  tau = zeros(J,N);
  for n = 1:N
    denom[n] = sum( exp( 
                -(alpha + piInc[K+1]*inc[n] + piAge[K+1]*age[n] + sigma[K+1]*v[n,K+1])*p
                + x[:,1:K]*(beta[1:K] + piInc[1:K]*inc[n] + piAge[1:K]*age[n] + Diagonal(sigma[1:K])*v'[1:K,n] )
                + xi )
                )                      
    tau[:,n] = exp(                           
                -(alpha + piInc[K+1]*inc[n] + piAge[K+1]*age[n] + sigma[K+1]*v[n,K+1])*p
                + x[:,1:K]*(beta[1:K] + piInc[1:K]*inc[n] + piAge[1:K]*age[n] + Diagonal(sigma[1:K])*v'[1:K,n] )
                + xi )/denom[n]
  end
  g[1:528] =  s - sum(tau[:,n], 2) / N
  g[529:538] = iv'param[31:558] - param[21:30]
end


function eval_jac_g(param, mode, rows, cols, values)

# Variables
  alpha = param[1];
  beta = param[2:5];
  piInc = param[6:10];
  piAge = param[11:15];
  sigma = param[16:20];
  g = param[21:30];
  xi = param[31:558];
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
                                          
    tau[:,n] = exp(                               
              -(alpha + piInc[K+1]*inc[n] + piAge[K+1]*age[n] + sigma[K+1]*v[n,K+1])*p
              + x[:,1:K]*(beta[1:K] + piInc[1:K]*inc[n] + piAge[1:K]*age[n] + Diagonal(sigma[1:K])*v'[1:K,n] )
              + xi )/denom[n]
  end

  # Define derivatives point by point
  for n = 1:N
    for j = 1:J
      d_alpha[j,n] = p[j]*tau[j,n] - sum(Diagonal(p)*tau[:,n])*tau[j,n]
      d_piInc[j,K+1,N] = (p[j]*tau[j,n] - sum(Diagonal(p)*tau[:,n])*tau[j,n])*inc[n]
      d_piAge[j,K+1,N] = (p[j]*tau[j,n] - sum(Diagonal(p)*tau[:,n])*tau[j,n])*age[n]
      d_sigma[j,K+1,N] = (p[j]*tau[j,n] - sum(Diagonal(p)*tau[:,n])*tau[j,n])*v'[K+1,n]
      for k=1:K
        d_beta[j,k,n] = ( x[j,k]*tau[j,n] - sum(Diagonal(x[:,k])*tau[:,n])*tau[j,n] )
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

  # Define derivative matrices
  D_alpha = (1/N)*sum(d_alpha, 2);
  D_beta = (1/N)*sum(d_beta, 3);
  D_piInc = (1/N)*sum(d_piInc, 3);
  D_piAge = (1/N)*sum(d_piAge, 3);
  D_sigma = (1/N)*sum(d_sigma, 3);
  D_theta = [D_beta D_alpha D_piInc D_piAge D_sigma];
  D_xi = (1/N)*sum(d_xi, 3);

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


#####################
##   Call Ipopt    ##
#####################

# Create the problem
prob = createProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess,
                     eval_f, eval_g, eval_grad_f, eval_jac_g)

# Set warm start and additional options
prob.x = zeros(n)
addOption(prob, "hessian_approximation", "limited-memory")
addOption(prob, "linear_solver", "ma57")

# Solve
status = solveProblem(prob)
println(Ipopt.ApplicationReturnStatus[status])
println(prob.x)
println(prob.obj_val)

