"""
- This script calculates the spinwave spectrum of the 120° antiferromagnetic order.
- Details of the formalism is summarized below.
"""
# ? Hamiltonian is given as
# ! H = ∑ₖ xₖ† Hₖ xₖ
# ? Nambu basis is selected as follows,
# ? Don't know the details, but sometimes they call it as "Nambu basis"
# ! xₖ = [ a₊ₖ  ]    xₖ† = [a₊ₖ†, a₋ₖ ]
# !      [ a₋ₖ† ]
# ? In this convention, hₖ is given as,
# ! hₖ = [ +A₊ₖ  -B₊ₖ ]
# !      [ -B₊ₖ  +A₋ₖ ] 
# ? compare the Hamiltonian with its original form.
# ? constant term in H₂ is omitted.
# ? for the beauty, we will use the following convention,
# ! H₂ = ∑ₖ S/4 [ Aₖ aₖ†aₖ + A₋ₖ a₋ₖ†a₋ₖ - Bₖ (a₋ₖaₖ + aₖ⁺a₋ₖ⁺) ]
# ? Bravais lattice case, we can easily calculate eigenstate energy of this problem.
# ! Ωₖ = ± √(AₖA₋ₖ - Bₖ²)


include("Param.jl");
using LinearAlgebra, GLMakie

# ! Be ware, in this example, 
# !  q vector is defined in terms of the reciprocal lattice vector.
# !  which means you don't need to multiply 2π to the q vector.
# !  It is different from different examples.

# as a baseline, we will start Colpa algorithm to 120° antiferromagnetic order
Qₘ = [ 1/3,  1/3];  
δQ₁ = [1E-6, 0E-6];
δQ₂ = [0E-6, 1E-6];

Qₘ₊₁ = Qₘ + δQ₁;  Qₘ₊₂ = Qₘ + δQ₂;
Qₘ₋₁ = Qₘ - δQ₁;  Qₘ₋₂ = Qₘ - δQ₂;

g = [[ 1  0]; 
     [ 0 -1]];
# ? commutation relation of xₖ and xₖ†

N = 100;
qList = [];  qAxis = [];


P₁ = [0.0,0.0];
P₂ = [1/2,0.0];
for x in range(0,stop=1,length=N)
  push!(qList, (1-x)*P₁+x*P₂)
  push!(qAxis, x)
end
qAxis = float.(qAxis)

Eₖ = zeros(N);  Ωₖ = zeros(N);

Qₛ = Qₘ₊₂;

for (it,q) in enumerate(qList)
  
  A₊ₖ = Aₖ(+q, J₁, J₂, J₃; Q₀ = Qₛ);
  A₋ₖ = Aₖ(-q, J₁, J₂, J₃; Q₀ = Qₛ);
  B₊ₖ = Bₖ(+q, J₁, J₂, J₃; Q₀ = Qₛ);
  
  h = [+A₊ₖ  -B₊ₖ; 
       -B₊ₖ  +A₋ₖ]

  # * Here, Colpa algorithm starts.
  # * Step 1: First perform Cholesky decomposition to hₖ
  
  h_safe = h + 1e-10*I;
    # ? For the numerical safety, 
    # ? we will add a small number to the diagonal elements of hₖ.
  
  C = cholesky(h_safe)
    # ? First, we perform Cholesky decomposition to hₖ.
    # ? hₖ = C.L * C.U
  
  # abs.((C.L*C.U) - h_safe) .< 10*eps(Float64)
    # ? You can check the result of the Cholesky decomposition.
    # ? It can be ckeced like this.
  
  UgL = C.U*g*C.L
    # ? Bogoliubov transformed basis version Hamiltonian is,
    # ? ⟹ U * g * L
    # ? Now you can use this matrix as normal number matrix.
    # ? Unitary transform solved the commutation problem.

  λs, Vs = eigen(UgL)
    # ? λs contains eigenvalues of UgL.
    # ? Vs contains eigenvectors of UgL.
    # ? UgL * Vs[k] = λs[k] * Vs[k]
  
  L = Vs' * UgL * Vs
    # ? This one is diagonal matrix generated from UgL.
    # ? Lᵢᵢ

  e = maximum(diag(L));

  Eₖ[it] = e;
  if A₊ₖ*A₋ₖ - B₊ₖ^2 > 0
    Ωₖ[it] = √(A₊ₖ*A₋ₖ - B₊ₖ^2);
  else 
    Ωₖ[it] = 0;    
  end
    # This one is a analytically obtained solution form.

end

f = Figure()
Axis(f[1, 1])
lines!(qAxis,Eₖ, color = :blue, linewidth = 4, linestyle = :solid, label = "Eₖ")
lines!(qAxis,Ωₖ, color = :red,  linewidth = 3, linestyle = :dash,  label = "Ωₖ")
f
