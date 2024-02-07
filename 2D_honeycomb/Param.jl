using LinearAlgebra

J₁ = +1.00;
dJ₁AtoB = [[ 0, 0], [+1, 0], [ 0,-1]];
dJ₁BtoA = [[ 0, 0], [-1, 0], [ 0,+1]];
S  = 1.0;

"""
It calculates the Fourier Transform of the exchange interaction.
We use the following convention:
  J(r) = ∑ₖ Jₖ exp(ikr)
"""
function Jₖ(k, J₁; dJ₁ = dJ₁AtoB)
  val = 0;
  for dᵢ in dJ₁
    val += J₁ * exp( 2*π*im*dot(k,dᵢ) );
  end
  return val;
end

function Aₖ(k, J₁, dJ₁; Q₀ = [0, 0])
  return 2*Jₖ(k, J₁; dJ₁=dJ₁) + 
           Jₖ(k+Q₀, J₁; dJ₁=dJ₁) + 
           Jₖ(k-Q₀, J₁; dJ₁=dJ₁) - 
         4*Jₖ(Q₀, J₁; dJ₁=dJ₁);
end
  
function Bₖ(k, J₁, dJ₁; Q₀ = [0, 0])
  return 2*Jₖ(k, J₁; dJ₁=dJ₁) - 
           Jₖ(k+Q₀, J₁; dJ₁=dJ₁) - 
           Jₖ(k-Q₀, J₁; dJ₁=dJ₁);
end