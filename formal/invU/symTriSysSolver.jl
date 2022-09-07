using LinearAlgebra


include("invSymTriSS.jl")

function solveTxb(T, b)
    n = length(b)
    xhat, yhat = invSymT2SS(T, false)
    triuxytb = xhat .* normCoef(yhat, b, false)
    trilyxtb = yhat .* vcat([0], normCoef(xhat[1:n-1], b[1:n-1], true))
    triuxytb + trilyxtb;
end

