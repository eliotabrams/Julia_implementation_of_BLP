# 2016 Winter Advanced IO PS2 
# Hyunmin Park, Eliot Abrams, Alexandre Sollaci

#=
julia_implementation_of_blp.jl

Julia code for implementing a BLP model using MPEC to solve for parameters
=#

#####################
##      Setup      ##
#####################

Pkg.add("Ipopt")
Pkg.add("JuMP")
using Ipopt
using JuMP
cd("/Users/eliotabrams/Desktop/Advanced\ Industrial\ Organization\ 2/Julia_implementation_of_BLP/")


#####################
##      Data       ##
#####################

# Load data
product = readdlm("dataset_cleaned.csv", ',');
product = convert(
    Array{Float64,2},product[2:size(product,1),1:size(product,2)]
    );

population = readdlm("population_data.csv", ',');
population = convert(
    Array{Float64,2},population[2:size(population,1),1:size(population,2)]
    );

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
tic()

# Setup the simple logit model
m = Model(solver = IpoptSolver(tol = 1e-8, max_iter = 1000, output_file = "logit.txt"));

# Define variables
@defVar(m, g[1:L]);
@defVar(m, xi[1:J]);
@defVar(m, alpha);
@defVar(m, beta[1:K]);

# We minimize the gmm objective with the identity as the weighting matrix
# subject to the constraints g = sum_j xi_j iv_j and market share equations
@setObjective(m, Min, sum{g[l]^2,l=1:L});
@addConstraint(
    m, 
    constr[l=1:L], 
    g[l]==sum{xi[j]*iv[j,l],j=1:J}
);
@addNLConstraint(
    m, 
    constr[j=1:J], 
    xi[j]==log(s[j])-log(s0[j])+alpha*p[j]-sum{beta[k]*x[j,k],k=1:K}
);

# Solve the model
status = solve(m);
print(status);
println("alpha = ", getValue(alpha));
println("beta = ", getValue(beta[1:K]));
toc()

# Save results to use in the setup of BLP Model
g_logit=getValue(g);
xi_logit=getValue(xi);
alpha_logit=getValue(alpha);
beta_logit=getValue(beta);


##########################
##      BLP Model       ##
##########################
tic()

# Calculate the optimal weighting matrix
W = inv((1/J)*iv'*Diagonal(diag(xi_logit*xi_logit'))*iv);

# Setup the BLP model
m = Model(solver = IpoptSolver(tol = 1e-8, max_iter = 1000, output_file = "BLP.txt"));

# Defining variables - set initial values to estimates from the logit model
# LATER WILL WANT TO RUN FOR MULTIPLE STARTING VALUES
@defVar(m, g[x=1:L], start=(g_logit[x]));
@defVar(m, xi[x=1:J], start=(xi_logit[x]));
@defVar(m, alpha, start=alpha_logit);
@defVar(m, beta[x=1:K], start=beta_logit[x]);

# Defining variables - heterogeneity parameters
@defVar(m, piInc[1:K+1]);
@defVar(m, piAge[1:K+1]);
@defVar(m, sigma[1:K+1]);

# We minimize the gmm objective - using the optimal weighting matrix! 
# subject to g = sum_j xi_j iv_j and market share equations - 
# Note that where we assign each shock could have minor effect on estimation results
# shock 1 : taste shock to constant
# shock 2 : taste shock to x1
# shock 3 : taste shock to x2
# shock 4 : taste shock to x3
# shock 5 : taste shock to price
@setObjective(m,Min,sum{sum{W[i,j]*g[i]*g[j],i=1:L},j=1:L});
@addConstraint(
    m, 
    constr[l=1:L], 
    g[l]==sum{xi[j]*iv[j,l],j=1:J}
);
@defNLExpr(
    denom[n=1:N],
    sum{exp(sum{(beta[k]+piInc[k]*inc[n]+piAge[k]*age[n]+sigma[k]*v[n,k])*x[h,k],k=1:K}
    -(alpha+piInc[K+1]*inc[n]+piAge[K+1]*age[n]+sigma[K+1]*v[n,K+1])*p[h]+xi[h]),h=1:J}
);
@addNLConstraint(
    m,
    constr[j=1:J], 
    s[j]==(1/N)*
        sum{exp(sum{(beta[k]+piInc[k]*inc[n]+piAge[k]*age[n]+sigma[k]*v[n,k])*x[j,k],k=1:K}
        -(alpha+piInc[K+1]*inc[n]+piAge[K+1]*age[n]+sigma[K+1]*v[n,K+1])*p[j]+xi[j])/denom[n]
       , n=1:N}
);


# Solve the model
status = solve(m)
print(status);
println("alpha = ", getValue(alpha))
println("beta = ", getValue(beta[1:K]))
println("piInc = ", getValue(piInc[1:K+1])
println("piAge = ", getValue(piAge[1:K+1])
println("sigma = ", getValue(sigma[1:K+1])

toc()

