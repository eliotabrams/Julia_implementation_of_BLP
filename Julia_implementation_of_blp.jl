"""
julia_implementation_of_blp.jl

Julia code for implementing a BLP model using MPEC to solve for parameters
"""

#####################
##      Setup      ##
#####################

Pkg.add("Optim")
using Optim


#####################
##   Simple Logit  ##
#####################

# Construct eval_f
function eval_f(param) 
    delta_j = log(data[:,2]) - log(data[:,14]) -  *(data[:,3:7] , param) 
    z = [data[:,3:6] data[:,8:13]]
    vector = *(delta_j', z)
    return *( *(vector, eye(10)), vector' )[1]
end 

# Read data
# id  share   const   x1  x2  x3  price   z1  z2  z3  z4  z5  z6  s0
# 1   2       3       4   5   6   7       8   9   10  11  12  13  14
data = readdlm("/Users/eliotabrams/Desktop/Advanced\ Industrial\ Organization\ 2/Julia_implementation_of_BLP/dataset_cleaned.csv", ',')
data = data[2:529,1:14]
data = convert(Array{Float64,2},data)
eval_f([1,1,1,1,1])

# Run GMM
results = optimize(eval_f, [0.0, 0.0, 0.0, 0.0])
results.minimum


prob = createProblem(
    5, 
    [-100.0, -100.0,-100.0,-100.0,-100.0], 
    [100.0, 100.0, 100.0, 100.0, 100.0], 
    1, 
    [-1.0], 
    [1.0], 
    0, 
    0,
    eval_f, 
    eval_g, 
    eval_grad_f, 
    0,
    )
addOption(prob, "hessian_approximation", "limited-memory")
1prob.x = [1.0, 5.0, 5.0, 1.0]
status = solveProblem(prob)

# ISSUE IS THAT WE NEED TO CALCULATE 
# 1) THE NUMBER OF NON-ZEROS,
# 2) GRADIENT OF THE CONSTRAINT
# 3) THE JACOBIAN

min gTWg
g

st. g = ts
st. shares = complex(function)

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





