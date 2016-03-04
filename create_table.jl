###############################
####      Assignment 3      ###
###############################
# Advanced IO 2
# Eliot Abrams, Hyunmin Park, Alexandre Sollaci


###############################
####         Setup          ###
###############################

using DataFrames
using StatsBase


###############################
####     Produce Table      ###
###############################

# Read in bootstrap results
boots = readdlm(string(2, "_estimate.csv"),',', header = true)[1];
for i in 3:10
boots = [boots; readdlm(string(i, "_estimate.csv"),',', header = true)[1]];
end
df = DataFrame(convert(Array{Float64}, boots[:,1:9]));
df[:x3] = round(df[:x3] ./ 100);

# Add standard errors to main results
results = readdlm(string(1, "_estimate.csv"),',', header = true)[1];
results = convert(Array{Float64}, results[:,1:9]);
for i in 4:9
    results = hcat(results, by(df, [:x3, :x2], df -> std(df[i]))[3]);
end
results = round(results, 3);

# Print results table
csvfile = open("results.csv","w");
write(csvfile, "Beta,Num_States,MPEC_theta1,MPEC_theta2,MPEC_RC,HM_theta1,HM_theta2,HM_RC");
for i in 1:size(results,2)
    write(csvfile, string("\n", results[i,2], ",", results[i,3], ","));
    for j in 4:9
        write(csvfile, string(results[i,j], " (",results[i,j+6],"),"));
    end
end
close(csvfile);


