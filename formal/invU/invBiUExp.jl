using LinearAlgebra

include("expNum.jl")
include("genTestMat.jl")

function expBidiagBackSubEn(U)
    n = size(U)[2]
    sx, xE = zeros(n), zeros(n)
    (sx[n], xE[n]) = expInv(fl2exp(U[n, n]))
    for i in n-1 : -1 : 1
        (sx[i], xE[i]) = expDivide(expTimes(fl2exp(-U[i, i+1]), (sx[i+1], xE[i+1])), fl2exp(U[i, i]))
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

