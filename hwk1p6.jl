
# Hwk1 Ps6 by Hyunmin Park
# Based on Newton's method on Rosenbrock's Function by Mihai on Jan 10, 2016


Pkg.add("JuMP")
Pkg.add("PyPlot")   

using JuMP
N=4
# N=10
# N=40
Iter=100
m=Model()
@defVar(m,x[1:N])
@setNLObjective(m,Min,sum{(x[i+1]-x[i])^2,i=1:N-1}+sum{(exp(x[i])-1)^2,i=1:N})

# define the abstract evaluator
d=JuMP.JuMPNLPEvaluator(m,JuMP.prepConstrMatrix(m))
println(d)

#  This portion demonstrates how to access the data from a jump model to include in the solver
using MathProgBase
import MathProgBase.MathProgSolverInterface
x=10*ones(N)
# initialize the Sparse Graph Structure
MathProgSolverInterface.initialize(d,[:Hess])
I,J=MathProgSolverInterface.hesslag_structure(d)
V=zeros(length(I))
g=zeros(N)
outGrad=zeros(Iter)
outGrad2=zeros(Iter)
# Newton Loop
for i=1:Iter
  MathProgSolverInterface.eval_grad_f(d,g,x)
  # print("Current gradient is ",g,"\n")
  # println("Norm of gradient is ",norm(g))
  # println(g)
  # println(x)
  outGrad[i]=log10(max(norm(g),1e-40))
    outGrad2[i]=norm(g)
  MathProgSolverInterface.eval_hesslag(d,V,x,1.0,Float64[])
  hess_raw=sparse(I,J,V)
  hess_raw=hess_raw+hess_raw'-spdiagm(diag(hess_raw))
  # above line: beccause eval_hesslag returns only 1/2 hessian
  # println(hess_raw)
  newtonDirection=-hess_raw\g
  # println(newtonDirection)
  x=x+newtonDirection
end
print(x)

using PyPlot
#PyPlot.pygui(true)
xgr = linspace(1, Iter-1, Iter-1)
#Quotient of adjacent gradient norms
gradQ = outGrad2[2:Iter]./outGrad2[1:Iter-1]
plot(xgr, gradQ, "b-", linewidth=2)
show()

using PyPlot
#PyPlot.pygui(true)
xgr = linspace(1, Iter-1, Iter-1)
#Quotient of adjacent gradient norms
gradQten = outGrad2ten[2:Iter]./outGrad2ten[1:Iter-1]
plot(xgr, gradQten, "b-", linewidth=2)
show()

using PyPlot
#PyPlot.pygui(true)
xgr = linspace(1, Iter-1, Iter-1)
#Quotient of adjacent gradient norms
gradQforty = outGrad2forty[2:Iter]./outGrad2forty[1:Iter-1]
plot(xgr, gradQforty, "b-", linewidth=2)
show()


