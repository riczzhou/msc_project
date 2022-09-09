function invSymT2SS(T, exactInv=true)
    C = cholesky(Matrix(T))
    U = Bidiagonal(C.U)
    (x, y) = invBidiagU(U,false)
    (xhat, yhat) = onepairSymSS(x, y)
    if exactInv
        S = generateSymSS(xhat, yhat)
        return S;
    end
    xhat, yhat;
end

function invSymT2SS_U(U, exactInv=true)
    (x, y) = invBidiagU(U,false)
    (xhat, yhat) = onepairSymSS(x, y)
    if exactInv
        S = generateSymSS(xhat, yhat)
        return S;
    end
    xhat, yhat;
end

function solverTxb(T, b)
    n = length(b)
    xhat, yhat = invSymT2SS(T, false)
    triuxytb = xhat .* normCoef(yhat, b, false)
    trilyxtb = yhat .* vcat([0], normCoef(xhat[1:n-1], b[1:n-1], true))
    triuxytb + trilyxtb;
end

function solver4seqTxb(xhat, yhat, b)
    n = length(b)
    triuxytb = xhat .* normCoef(yhat, b, false)
    trilyxtb = yhat .* vcat([0], normCoef(xhat[1:n-1], b[1:n-1], true))
    triuxytb + trilyxtb;
end