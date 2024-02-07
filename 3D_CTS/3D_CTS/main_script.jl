using GLMakie
using LinearAlgebra
using Statistics
include("Func_bundle.jl");

println("====================================================================================")
println("This code will calculate the optimal magnetic structure for the given parameters.")
println("Since we are dealing with ABAB stacking of triangular lattice,")
println("optimal magnetic structure should be specified by,")
println("the wave vector Qₘ and the angle θ between the two sublattices.")

qList = range(0,1/3,length=100);
Es = zeros(size(qList)[1],3);
vs = zeros(Complex{Float64},size(qList)[1],2,2);

for (it,q) in enumerate(qList)
  
  Q = [q 0 0];
  
  Jmat = zeros(Complex{Float64},2,2);
  
  Jmat[1,1] += Jₖ(Q,J₁,dJ₁AtoA);    Jmat[2,2] += Jₖ(Q,J₁,dJ₁BtoB);
  Jmat[1,1] += Jₖ(Q,J₂,dJ₂AtoA);    Jmat[2,2] += Jₖ(Q,J₂,dJ₂BtoB);
  Jmat[1,1] += Jₖ(Q,J₃,dJ₃AtoA);    Jmat[2,2] += Jₖ(Q,J₃,dJ₃BtoB);
  Jmat[1,2] += Jₖ(Q,Jc₁,dJc₁AtoB);  Jmat[2,1] += Jₖ(Q,Jc₁,dJc₁BtoA);
  Jmat[1,2] += Jₖ(Q,Jc₂,dJc₂AtoB);  Jmat[2,1] += Jₖ(Q,Jc₂,dJc₂BtoA);
  
  λs,Vs = eigen(Jmat)
  Es[it,1:2] = real(λs);
  Es[it,3] = real(mean(diag(Jmat)));
  vs[it,:,:] = Vs;
end

idx = argmin(Es[:,1])
qOpt = qList[idx];
vecTmp = vs[idx,:,1];
vec = [conj(vecTmp[1])/norm(conj(vecTmp[1]))*vecTmp[1], 
       conj(vecTmp[1])/norm(conj(vecTmp[1]))*vecTmp[2]]
deg = angle(vec[2])/2π*360
qOptApprox = round(qOpt,digits=4)
degApprox = round(deg,digits=2)

println("====================================================================================")
println("Optimal ordering vector is Q = $(qOptApprox) aˣ")
println("Optimal relative angle b/w A and B sublattice is θ ~ $(degApprox) degree")


fig = Figure(size = (800, 600));
ax = Axis(fig[1,1], xlabel = "q", ylabel = "E(q)",xticks = (0:1/9:1/3, ["0", "1/9", "2/9", "1/3"]));
lines!(ax, qList,Es[:,1],color=:blue);
lines!(ax, qList,Es[:,2],color=:red);
lines!(ax, qList,Es[:,3],color=:green);
fig