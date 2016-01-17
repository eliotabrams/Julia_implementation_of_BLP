"""
optimization_call.jl

Code for implementing a BLP model using MPEC to solve for parameters
"""

#####################
## Import Packages ##
#####################
using ipopt
using Calculus
Pkg.add("DataFrames")
using DataFrames
using Ipopt


#####################
##  Code Snippets  ##
#####################

# Set working directory
cd("/Users/eliotabrams/Desktop/Advanced\ Industrial\ Organization\ 2/Problem\ Set\ 2/")


# Printing
println("hello world")

# Function declaration
function sphere_vol(r)
    # julia allows Unicode names (in UTF-8 encoding)
    # so either "pi" or the symbol π can be used
    return 4/3*pi*r^3
end

# Function calling
sphere_vol(10)

# Data importation
basic_table = readtable("dataset_cleaned.csv",separator=',',header=true)

# Data work
basic_table[1, :price]
basic_table[basic_table[:price] .> 1, :]
basic_table[:new_column] = map(x -> x - 10, basic_table[:price])


