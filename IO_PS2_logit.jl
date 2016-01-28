# 2016 Winter Advanced IO PS2 by Eliot Abrams, Hyunmin Park, Alexandre Sollaci

tic()

# Load data
using DataFrames
product = readtable("dataset_cleaned.csv")
population = readtable("population_data.csv")

# Define variables
x = product[:,3:6]
p = product[:,7]
z = product[:,8:13]
s0 = product[:,14]
s = product[:,2]
iv = hcat(x,z)

# Setting up the model
using JuMP
K = size(x,2)
L = K+size(z,2)
J = size(x,1)
m = Model()


# Defining variables
@defVar(m, g[1:L])
@defVar(m, xi[1:J])
@defVar(m, alpha)
@defVar(m, beta[1:K])
@setObjective(m,Min,sum{g[l]^2,l=1:L})

# g = sum_j xi_j iv_j
l = 1
while l <= L
	@addConstraint(m,g[l]==sum{xi[j]*iv[j,l],j=1:J}) 
	l += 1
end

# market share equations (loop)
j = 1
while j <= J
	@addNLConstraint(m,xi[j]==log(s[j])-log(s0[j])+alpha*p[j]-sum{beta[k]*x[j,k],k=1:K})
	j += 1
end

using Ipopt
setSolver(m,IpoptSolver(tol = 1e-10, max_iter = 200, output_file = "results.txt"))

#for l=1:L
#  setValue(g[l],1)
#end

status = solve(m)
print(status)

println("alpha = ", getValue(alpha))
println("beta = ", getValue(beta[1:K]))

toc()
