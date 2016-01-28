
#####################
##      Setup      ##
#####################

Pkg.add("Ipopt")
using Ipopt


#####################
##   Simple Logit  ##
#####################
tic()

# Construct eval_f
function eval_f(param) 
    delta_j = log(data[:,2]) - log(data[:,14]) -  *(data[:,3:7] , param) 
    z = [data[:,3:6] data[:,8:13]]
    vector = *(delta_j', z)
    return *( *(vector, eye(10)), vector' )[1]
end 

# Construct eval_g
function eval_g(param, g)
    return g[1] = 0
end

# Construct eval_grad_f NEEDS TO BE FIXED!!!
# DANG ITS HARD TO CONSTRUCT THE GRADIENTS BY HAND :-(
function eval_grad_f(param, grad_f)
  grad_f[1] = *(param', *(*(data[:,3:7]', [data[:,3:6] data[:,8:13]]),*(data[:,3:7]', [data[:,3:6] data[:,8:13]])'))[1]
  grad_f[2] = *(param', *(*(data[:,3:7]', [data[:,3:6] data[:,8:13]]),*(data[:,3:7]', [data[:,3:6] data[:,8:13]])'))[2]
  grad_f[3] = *(param', *(*(data[:,3:7]', [data[:,3:6] data[:,8:13]]),*(data[:,3:7]', [data[:,3:6] data[:,8:13]])'))[3]
  grad_f[4] = *(param', *(*(data[:,3:7]', [data[:,3:6] data[:,8:13]]),*(data[:,3:7]', [data[:,3:6] data[:,8:13]])'))[4]
  grad_f[5] = *(param', *(*(data[:,3:7]', [data[:,3:6] data[:,8:13]]),*(data[:,3:7]', [data[:,3:6] data[:,8:13]])'))[5]
end

# Construct sparsity structure of the Jacobian
function eval_jac_g(param, mode, rows, cols, values)
  if mode == :Structure
    rows[1] = 1; cols[1] = 1
  end
end

# Read data
# id  share   const   x1  x2  x3  price   z1  z2  z3  z4  z5  z6  s0
# 1   2       3       4   5   6   7       8   9   10  11  12  13  14
data = readdlm("/Users/eliotabrams/Desktop/Advanced\ Industrial\ Organization\ 2/Julia_implementation_of_BLP/dataset_cleaned.csv", ',')
data = data[2:529,1:14]
data = convert(Array{Float64,2},data)
# eval_f([1,1,1,1,1])

# Run GMM
prob = createProblem(
    5, 
    [-100.0, -100.0,-100.0,-100.0,-100.0], 
    [100.0, 100.0, 100.0, 100.0, 100.0], 
    1, 
    [-0.0], 
    [0.0], 
    0, 
    0,
    eval_f, 
    eval_g, 
    eval_grad_f, 
    eval_jac_g,
    )
addOption(prob, "hessian_approximation", "limited-memory")
prob.x = [1.0, 1.0, 5.0, 5.0, 1.0]
status = solveProblem(prob)
println(Ipopt.ApplicationReturnStatus[status])
println(prob.x)
println(prob.obj_val)

toc()



