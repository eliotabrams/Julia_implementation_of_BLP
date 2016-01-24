"""
julia_implementation_of_blp.jl

Julia code for implementing a BLP model using MPEC to solve for parameters
"""

#####################
##      Setup      ##
#####################

Pkg.add("Ipopt")
Pkg.add("Optim")

using Ipopt
using Optim

cd("/Users/eliotabrams/Desktop/Advanced\ Industrial\ Organization\ 2/Julia_implementation_of_BLP/")


#####################
##   Simple Logit  ##
#####################

# Construct f
function f(param, share, s0, char, instr)
    return (log(share) - log(s0) -  dot(char, param)) * cat(1, char[1:3], instr)
end

# Construct objective 
function objective(param)
    g = zeros(9)
    for i = 2:size(data)[1]
        g += f(
                param,
                data[i,:][2], 
                data[i,:][14], 
                data[i,:][4:7], 
                data[i,:][8:13]
            )
    end
    return dot(g, g) 
end

# Run GMM
data = readdlm("dataset_cleaned.csv", ',');
objective([1,1,1,1])
results = optimize(objective, [0.0, 0.0, 0.0, 0.0])
results.minimum


# Import data
# id  share   const   x1  x2  x3  price   z1  z2  z3  z4  z5  z6  s0
# 1   2       3       4   5   6   7       8   9   10  11  12  13  14
#share = .2
#s0 = .1
#char = [2, 3, 4, 2]
#param = [1, 1, 1, 2]
#instr = [1, 2, 3, 4, 5, 6]
#f(param, share, s0, char, instr)
#objective(param)

