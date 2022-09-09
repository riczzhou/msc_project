function invSymT2SSexp(T, exactInv=true)
    C = cholesky(Matrix(T))
    U = Bidiagonal(C.U)
    (sx, xE), (sy, yE) = invBiUexp(U,false)
    (sxh, xhE), (syh, yhE) = onepairSymSSexp((sx, xE), (sy, yE))
    if exactInv
        S = generateSymSSexp((sxh, xhE), (syh, yhE))
        return S;
    end
    (sxh, xhE), (syh, yhE);
end

function invSymT2SSexp_U(U, exactInv=true)
    (sx, xE), (sy, yE) = invBiUexp(U,false)
    (sxh, xhE), (syh, yhE) = onepairSymSSexp((sx, xE), (sy, yE))
    if exactInv
        S = generateSymSSexp((sxh, xhE), (syh, yhE))
        return S;
    end
    (sxh, xhE), (syh, yhE);
end

function solverTxbexp(T, b)
    n = length(b)
    sb, bE = fl2exp(b)
    (sxh, xhE), (syh, yhE) = invSymT2SSexp(T, false)
    striuxytb, triuxytbE = expTimes((sxh, xhE), normCoefexp((syh, yhE), (sb, bE), false))
    strilyxtb_, trilyxtbE_ = zeros(n-1), zeros(n-1)
    sb_, bE_ = sb[1:n-1], bE[1:n-1]
    sxh_, xhE_ = sxh[1:n-1], xhE[1:n-1]
    syh_, yhE_ = syh[2:n], yhE[2:n]
    strilyxtb_, trilyxtbE_ = expTimes((syh_, yhE_), normCoefexp((sxh_, xhE_), (sb_, bE_), true))
    x = zeros(n)
    x[1] = exp2fl((striuxytb[1], triuxytbE[1]))
    x[2:n] = exp2fl((striuxytb[2:n], triuxytbE[2:n])) + exp2fl((strilyxtb_, trilyxtbE_))
    x;
end

function solver4seqTxbexp((sxh, xhE), (syh, yhE), b)
    n = length(b)
    sb, bE = fl2exp(b)
    striuxytb, triuxytbE = expTimes((sxh, xhE), normCoefexp((syh, yhE), (sb, bE), false))
    strilyxtb_, trilyxtbE_ = zeros(n-1), zeros(n-1)
    sb_, bE_ = sb[1:n-1], bE[1:n-1]
    sxh_, xhE_ = sxh[1:n-1], xhE[1:n-1]
    syh_, yhE_ = syh[2:n], yhE[2:n]
    strilyxtb_, trilyxtbE_ = expTimes((syh_, yhE_), normCoefexp((sxh_, xhE_), (sb_, bE_), true))
    x = zeros(n)
    x[1] = exp2fl((striuxytb[1], triuxytbE[1]))
    x[2:n] = exp2fl((striuxytb[2:n], triuxytbE[2:n])) + exp2fl((strilyxtb_, trilyxtbE_))
    x;
end