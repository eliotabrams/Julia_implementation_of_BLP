# 2016 Winter Advanced IO PS2 
# Hyunmin Park, Eliot Abrams, Alexandre Sollaci

#=
julia_implementation_of_blp.jl

Julia code for implementing a BLP model using MPEC to solve for parameters
=#

#####################
##      Setup      ##
#####################

#Pkg.add("Ipopt")
#Pkg.add("JuMP")
using Ipopt
using JuMP
using DataFrames
#cd("/Users/eliotabrams/Desktop/Advanced\ Industrial\ Organization\ 2/Julia_implementation_of_BLP")
EnableNLPResolve()

#####################
##      Data       ##
#####################

# Load data
product = DataFrames.readtable("small_dataset_cleaned.csv", separator = ',', header = true);
population = DataFrames.readtable("small_population_data.csv", separator = ',', header = true);

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
print(status)
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
iv = convert(Array, iv)
W = inv((1/J)*iv'*Diagonal(diag(xi_logit*xi_logit'))*iv);

# Setup the BLP model
BLP = Model(solver = IpoptSolver(hessian_approximation = limited-memory, tol = 1e-6, max_iter = 1000, output_file = "BLP.txt"));

# Defining variables - set initial values to estimates from the logit model
@defVar(BLP, g[x=1:L], start=(g_logit[x]));
@defVar(BLP, xi[x=1:J], start=(xi_logit[x]));
@defVar(BLP, alpha, start=alpha_logit);
@defVar(BLP, beta[x=1:K], start=beta_logit[x]);

# Defining variables - heterogeneity parameters
@defVar(BLP, piInc[1:K]);
@defVar(BLP, piAge[1:K]);
@defVar(BLP, sigma[1:K]);

# We minimize the gmm objective - using the optimal weighting matrix! 
# subject to g = sum_j xi_j iv_j and market share equations - 
# Note that where we assign each shock could have minor effect on estimation results
# shock 1 : taste shock to price
# shock 2 : taste shock to x1
# shock 3 : taste shock to x2
# shock 4 : taste shock to x3
@setObjective(BLP,Min,sum{sum{W[i,j]*g[i]*g[j],i=1:L},j=1:L});
@addConstraint(
    BLP, 
    constr[l=1:L], 
    g[l]==sum{xi[j]*iv[j,l],j=1:J}
);
@defNLExpr(
    BLP,
    denom[n=1:N],
    sum{
        exp(beta[1]
            -(alpha+piInc[1]*inc[n]+piAge[1]*age[n]+sigma[1]*v[n,1])*p[h]
            +sum{(beta[k]+piInc[k]*inc[n]+piAge[k]*age[n]+sigma[k]*v[n,k])*x[h,k],k=2:K}
            +xi[h]
            )
    , h=1:J}
);
@addNLConstraint(
    BLP,
    constr[j=1:J], 
    s[j]==(1/N)*
        sum{
            exp(beta[1]
                -(alpha+piInc[1]*inc[n]+piAge[1]*age[n]+sigma[1]*v[n,1])*p[j]
                +sum{(beta[k]+piInc[k]*inc[n]+piAge[k]*age[n]+sigma[k]*v[n,k])*x[j,k],k=2:K}
                +xi[j]
            )/denom[n]
       ,n=1:N}
);


status = solve(BLP);

# Print the results
print(status)
println("alpha = ", getValue(alpha))
println("beta = ", getValue(beta[1:K]))
println("piInc = ", getValue(piInc[1:K])
println("piAge = ", getValue(piAge[1:K])
println("sigma = ", getValue(sigma[1:K])


