#=
2016 Winter Advanced IO PS2 
Hyunmin Park, Eliot Abrams, Alexandre Sollaci

Julia code for implementing a BLP model using MPEC to solve for parameters
=#

#####################
##      Setup      ##
#####################

using Ipopt
using KNITRO
using JuMP
using DataFrames
EnableNLPResolve()


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
logit = Model(solver = IpoptSolver(tol = 1e-8, max_iter = 1000, output_file = "logit.txt"));

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
println("alpha = ", getValue(alpha))
println("beta = ", getValue(beta[1:K]))

# Save results to use in the setup of BLP Model
g_logit=getValue(g);
xi_logit=getValue(xi);
alpha_logit=getValue(alpha);
beta_logit=getValue(beta);


##########################
##      BLP Model       ##
##########################

# Calculate the optimal weighting matrix
iv = convert(Array, iv);
W = inv((1/J)*iv'*Diagonal(diag(xi_logit*xi_logit'))*iv);

# Setup the BLP model
BLP = Model(solver = KnitroSolver(KTR_PARAM_HESSOPT=6, KTR_PARAM_OUTMODE=2, KTR_PARAM_LINSOLVER=5, KTR_PARAM_MAXIT=30));

# Defining variables - set initial values to estimates from the logit model
@defVar(BLP, g[x=1:L], start=(g_logit[x]));
@defVar(BLP, xi[x=1:J], start=(xi_logit[x]));
@defVar(BLP, alpha, start=alpha_logit);
@defVar(BLP, beta[x=1:K], start=beta_logit[x]);

# Defining variables - heterogeneity parameters
@defVar(BLP, piInc[1:K+1]);
@defVar(BLP, piAge[1:K+1]);
@defVar(BLP, sigma[1:K+1]);

# We minimize the gmm objective - using the optimal weighting matrix
# subject to g = sum_j xi_j iv_j and market share equations - 
# Note that where we assign each shock could have minor effect on estimation results
# shock 1 : taste shock to price
# shock 2 : taste shock to x1
# shock 3 : taste shock to x2
# shock 4 : taste shock to x3
# shock 5 : taste shock to constant
@setObjective(BLP,Min,sum{sum{W[i,j]*g[i]*g[j],i=1:L},j=1:L});
@addConstraint(
    BLP, 
    constr[l=1:L], 
    g[l] == sum{xi[j]*iv[j,l],j=1:J}
);

# Trick to increase the sparsity and aid AD
@defVar(BLP, summand[1:N,1:J])
@addConstraint(BLP,
              summand_constr[n=1:N,h=1:J],
              summand[n,h] == -(alpha+piInc[K+1]*inc[n]+piAge[K+1]*age[n]+sigma[K+1]*v[n,K+1])*p[h]
              +sum{(beta[k]+piInc[k]*inc[n]+piAge[k]*age[n]+sigma[k]*v[n,k])*x[h,k],k=1:K}
              +xi[h]
);
@defNLExpr(
    BLP,
    denom[n=1:N],
    sum{exp(summand[n,h]), h=1:J}
);
@addNLConstraint(
    BLP,
    constr[j=1:J], 
    s[j]==(1/N)*sum{exp(summand[n,j])/denom[n],n=1:N}
);

# Solve the model
status = solve(BLP);

# Print the results
println("alpha = ", getValue(alpha))
println("beta = ", getValue(beta[1:K]))
println("piInc = ", getValue(piInc[1:K+1]))
println("piAge = ", getValue(piAge[1:K+1]))
println("sigma = ", getValue(sigma[1:K+1]))



