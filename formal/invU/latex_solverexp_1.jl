function backSubVecExp(U, b)
    T = eltype(U)
    n = size(U)[2]
    sx, xE = zeros(T, n), zeros(T, n)
    (sx[n], xE[n]) = expDivide(fl2exp(b[n]), fl2exp(U[n, n]))

    for i in n-1 : -1 : 1
        sum_ = zero(T)
        for j in i+1 : n
            sum_ += exp2fl(expTimes(fl2exp(U[i, j]), (sx[j], xE[j])))
        end
        (sx[i], xE[i]) = expDivide(fl2exp(b[i] - sum_), fl2exp(U[i, i]))
    end
    sx, xE;
end

function bandedBackSubVecExp(U, b, bw)
    T = eltype(U)
    n = size(U)[2]
    sx, xE = zeros(T, n), zeros(T, n)
    (sx[n], xE[n]) = expDivide(fl2exp(b[n]), fl2exp(U[n, n]))
    for i in n-1 : -1 : 1
        sum_ = zero(T)
        for j in i+1 : min(n, i+bw)
            sum_ += exp2fl(expTimes(fl2exp(U[i, j]), (sx[j], xE[j])))
        end
        (sx[i], xE[i]) = expDivide(fl2exp(b[i] - sum_), fl2exp(U[i, i]))
    end
    sx, xE;
end

function invBiUexp(U, exactInv=true)
    n = size(U)[2]
    (sx, xE) = expBidiagBackSubEn(U)
    (sy, yE) = expInv(expTimes(fl2exp(U[diagind(U)]), (sx, xE)))
    if exactInv
        Uinv = zeros(n, n)
        for i in 1:n
            for j in i:n
                Uinvsi, UinvEi = expTimes((sx[i], xE[i]), (sy[j], yE[j]))
                Uinv[i, j] = exp2fl((Uinvsi, UinvEi))
            end
        end
        return Uinv;
    end
    (sx, xE), (sy, yE);
end