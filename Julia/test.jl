using DataFrames
using JuMP
using Gurobi




m = Model(solver=GurobiSolver())

@variable(m,x[1:3])
@objective(m,Min,x'*4*eye(3)*x+x'*4*eye(3)*x)
@constraint(m,x.>=0)
@constraint(m,sum(x)>=10)

status = solve(m)
println(getvalue(x))

break

m = Model(solver=GurobiSolver())

@variable(m,x1)
@variable(m,x2)
@objective(m,Min,x1^2+4*x2^2)
@constraint(m,x1>=0)
@constraint(m,x2>=0)
@constraint(m,x2+x1>=10)

status = solve(m)
