"""
julia_implementation_of_blp.jl

Julia code for implementing a BLP model using MPEC to solve for parameters
"""

#####################
##      Setup      ##
#####################

Pkg.add("Ipopt")
Pkg.add("JuMP")
Pkg.add("ReverseDiffSource")
Pkg.add("Optim")

using Ipopt
using JuMP
using ReverseDiffSource
using Optim

cd("/Users/eliotabrams/Desktop/Advanced\ Industrial\ Organization\ 2/Julia_implementation_of_BLP/")


#####################
##   Simple Logit  ##
#####################

# Construct f
function eval_f(param) 
    delta_j = log(data[:,2]) - log(data[:,14]) -  *(data[:,4:7] , param) 
    g =   [
        sum(delta_j .* data[:,4]),
        sum(delta_j .* data[:,5]),
        sum(delta_j .* data[:,6]),
        sum(delta_j .* data[:,8]),
        sum(delta_j .* data[:,9]),
        sum(delta_j .* data[:,10]),
        sum(delta_j .* data[:,11]),
        sum(delta_j .* data[:,12]),
        sum(delta_j .* data[:,13])
    ]
    return dot(g, g)
end 

# id  share   const   x1  x2  x3  price   z1  z2  z3  z4  z5  z6  s0
# 1   2       3       4   5   6   7       8   9   10  11  12  13  14
data = readdlm("dataset_cleaned.csv", ',')
data = data[2:529,1:14]
data = convert(Array{Float64,2},data)

eval_f([1,1,1,1])

results = optimize(eval_f, [0.0, 0.0, 0.0, 0.0])
results.minimum



# Run GMM
g = rdiff(eval_f, (ones(4),), order=1)
prob = createProblem(4, objective)
prob.x = [1.0, 5.0, 5.0, 1.0]
status = solveProblem(prob)
prob = createProblem(
    4, 
    [-100.0,-100.0,-100.0,-100.0], 
    [100.0, 100.0, 100.0, 100.0], 
    1, 
    [0.0], 
    [0.0], 
    0, 
    0,
    objective, 
    0, 
    eval_grad_f, 
    0,
    0
    )
addOption(prob, "hessian_approximation", "limited-memory")
prob = createProblem(n, x_L, x_U, m, g_L, g_U, 8, 10,
                     eval_f, eval_g, eval_grad_f)

# ISSUE IS THAT WE NEED TO CALCULATE 
# 1) THE NUMBER OF NON-ZEROS,
# 2) GRADIENT OF THE CONSTRAINT
# 3) THE JACOBIAN
# 4) THE SPARSITY MATRICES

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





