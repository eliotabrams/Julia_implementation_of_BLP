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

using Ipopt
using JuMP
using ReverseDiffSource

cd("/Users/eliotabrams/Desktop/Advanced\ Industrial\ Organization\ 2/Julia_implementation_of_BLP/")


#####################
##   Simple Logit  ##
#####################

# Construct f
function f(param, share, s0, char, instr)
    return log(share) - log(s0) -  *(char, param)  [char instr]
end


log(share) - log(s0) -  *(char, param)  [char instr]


# Construct objective 
function objective(param)
    g = zeros(10)
    for i = 2:size(data)[1]
        g += f(
                param,
                data[i,:][2], 
                data[i,:][14], 
                data[i,:][4:7], 
                data[i,:][8:13]
            )
    end
    return sum(g.*g) 
end


param
share = data[:,2]
s0 = data[:,14]
char = data[:,4:7] 
instr = data[:,8:13]


# Import data
# id  share   const   x1  x2  x3  price   z1  z2  z3  z4  z5  z6  s0
# 1   2       3       4   5   6   7       8   9   10  11  12  13  14
data = readdlm("dataset_cleaned.csv", ',')
data = data[2:529,1:14]
data = convert(Array{Float64,2},data)


share = .2
s0 = .1
char = [2, 3, 4, 2]
param = [1, 1, 1, 1]
instr = [1, 2, 3, 4, 5, 6]
f(param, share, s0, char, instr)
objective(param)

# Run GMM

using DualNumbers

g = rdiff(objective, (ones(4),), order=1)

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





