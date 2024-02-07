include("Param.jl");
using GLMakie

Qₘ = [0,0];

qList = [];

q₁sweep = range(-1,+1,length=80);
q₂sweep = range(-1,+1,length=80);

for q₁ in q₁sweep, q₂ in q₂sweep
  push!(qList,[q₁,q₂]);
end

# for q in range(0,1/2,length=100)
#   push!(qList,[0,q]);
# end

Ωₖ = zeros(size(qList)[1],2);
qList = float.(qList);


g = [+1  0  0  0;
      0 +1  0  0;
      0  0 -1  0;
      0  0  0 -1;];

# ! Colpa method starts from here.
# ! in this case, honeycomb case, we only have J₁ antiferromagnet
# ! so we can simplify the code a bit
# ! anyway, Neel state is ground state.

for (idx,q) in enumerate(qList)

  hₖ = zeros(Complex{Float64},4,4);
  for i in 1:4
    hₖ[i,i] = 3*J₁*S;
  end

  hₖ[1,4] = -S * Jₖ(+q, J₁; dJ₁ = dJ₁BtoA);
  hₖ[2,3] = -S * Jₖ(+q, J₁; dJ₁ = dJ₁AtoB);
  hₖ[3,2] = -S * Jₖ(-q, J₁; dJ₁ = dJ₁AtoB);
  hₖ[4,1] = -S * Jₖ(-q, J₁; dJ₁ = dJ₁BtoA);

  # h_safe = 0.5*(hₖ+hₖ') + 1E-8I;
  h_safe = hₖ + 1E-8I;

  if ishermitian(h_safe) == false
    for row in eachrow(h_safe)
      println(row)
    end
    error("hₖ is not hermitian")
  end

  C = cholesky(h_safe);

  UgL = C.U*g*C.L;

  λs, Vs = eigen(UgL)

  L = Vs' * UgL * Vs
  # println(L)

  e = sort(λs);
  Ωₖ[idx,1] = e[4];  Ωₖ[idx,2] = e[3];
end

# using GLMakie.lines, plot the dispersion relation

basis = [+1.0 +0.5;
          0.0 √3/2 ]
x = [];
y = [];

for q in qList
  tmp = basis*q;
  push!(x,tmp[1]);  push!(y,tmp[2]);
end

x = float.(x);  y = float.(y);

f = Figure()
ax = Axis(f[1, 1],aspect=AxisAspect(√3))
scatter!(ax,x,y,color=Ωₖ[:,1]);
# xlims!(ax,-1.5,1/5)
# ylims!(ax,-√3/2,√3/2)
f