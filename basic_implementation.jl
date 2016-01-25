
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
results = optimize(eval_f, [0.0, 0.0, 0.0, 0.0, 0.0])
results.minimum

