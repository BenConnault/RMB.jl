using Base.Test
using RMB

# println("hi")


tol=1e-10

a=rand(2,5)
b=rand(5,2)
@test norm(I-RMB.explicitinv(a,b)*[a zeros(2,2);eye(5) b]) < tol 


