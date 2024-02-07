include("Param.jl")

function Jₖ(k, Jₓ, dJₓ)

  Jsum = 0;
  for d in dJₓ
    Jsum += Jₓ * exp(2π*im*dot(k,d));
  end

  return Jsum;
end