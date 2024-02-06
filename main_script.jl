include("Param.jl");
using LinearAlgebra
using GLMakie

# ! Be ware, 
# !  In this example, 
# !  q vector does not include the factor of 2π.
# !  It is just a simple vector in the reciprocal space.
# !  It is different from different examples.

# as a baseline, we will start Colpa algorithm to 120° antiferromagnetic order
Q₀ = [1/3,1/3];

# Hamiltonian is given as
# ! H = ∑ₖ xₖ† Hₖ xₖ

# Nambu basis is selected as follows,
# Don't know the details, but sometimes they call it as "Nambu basis"
# ! xₖ  = [a₊ₖ , a₋ₖ†]ᵀ
# ! xₖ† = [a₊ₖ†, a₋ₖ ]

# In this convention, hₖ is given as,
# ! hₖ = 
# !      [ +A₊ₖ  -B₊ₖ ]
# !      [ -B₊ₖ  +A₋ₖ ] 

# compare the Hamiltonian with its original form.
# constant term in H₂ is omitted.
# for the beauty, we will use the following convention,
# ? H₂ = ∑ₖ S/4 [ Aₖ aₖ†aₖ + A₋ₖ a₋ₖ†a₋ₖ - Bₖ (a₋ₖaₖ + aₖ⁺a₋ₖ⁺) ]

q = [1/7,0];

A₊ₖ = Aₖ(+q, J₁, J₂, J₃; Q₀ = Q₀);
A₋ₖ = Aₖ(-q, J₁, J₂, J₃; Q₀ = Q₀);
B₊ₖ = Bₖ(+q, J₁, J₂, J₃; Q₀ = Q₀);

h = [+A₊ₖ  -B₊ₖ; 
     -B₊ₖ  +A₋ₖ]

eigen!(h)

g = [1 0; 0 -1];
# ? commutation relation of xₖ and xₖ†

# * Here, Colpa algorithm starts.

# * Step 1: First perform Cholesky decomposition to hₖ

h_safe = h + 1e-10*I;
# ? For the numerical safety, we will add a small number to the diagonal elements of hₖ.

C = cholesky(h_safe)

abs.((C.L*C.U) - h_safe) .< 10*eps(Float64)
# ? cholesky decomposition is successful.

g_sandwiched = C.U*g*C.L

Λ, U = eigen(g_sandwiched)

U_tmp = zeros(size(U))
sz = size(U);
for (idx,u) in enumerate(eachcol(U))
  U_tmp[:,sz[2]-idx+1] = u
end
U = U_tmp;

L = U' * g_sandwiched * U
