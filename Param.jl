using LinearAlgebra

J₁ = +1.00;
dJ₁ = [[+1, 0], [ 0,+1], [+1,+1]];
J₂ = +0.00;
dJ₂ = [[+2,+1], [+1,+2], [+1,-1]];
J₃ = +0.00;
dJ₃ = [[+2, 0], [ 0,+2], [+2,+2]];
S  = 1.5;

"""
It calculates the Fourier Transform of the exchange interaction.
We use the following convention:
  J(r) = ∑ₖ Jₖ exp(ikr)
"""
function Jₖ(k, J₁, J₂, J₃; dJ₁ = dJ₁, dJ₂ = dJ₂, dJ₃ = dJ₃)
  val = 0;
  for (J,d) in zip([J₁,J₂,J₃],[dJ₁,dJ₂,dJ₃])
    for dᵢ in d
      val += J * cos( 2*π*dot(k,dᵢ) );
    end
  end
  return val;
end

function Aₖ(k, J₁, J₂, J₃; Q₀ = [acos(-1/2 * (1+J₂)/(J₂+2*J₃)), 0])
  return 2*Jₖ(k, J₁, J₂, J₃) + Jₖ(k+Q₀, J₁, J₂, J₃) + Jₖ(k-Q₀, J₁, J₂, J₃) - 4*Jₖ(Q₀, J₁, J₂, J₃);
end
  
function Bₖ(k, J₁, J₂, J₃; Q₀ = [acos(-1/2 * (1+J₂)/(J₂+2*J₃)), 0])
  return 2*Jₖ(k, J₁, J₂, J₃) - Jₖ(k+Q₀, J₁, J₂, J₃) - Jₖ(k-Q₀, J₁, J₂, J₃);
end
  
