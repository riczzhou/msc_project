using LinearAlgebra

include("expNum.jl")
include("genTestMat.jl")

function expBidiagBackSubEn(U)
    n = size(U)[2]
    T = eltype(U)
    sx, xE = zeros(T, n), zeros(T, n)
    (sEn, EnE) = T.((1.0, 0.0))
    (sx[n], xE[n]) = expInv(fl2exp(U[n, n]))
    for i in n-1 : -1 : 1
        (sx[i], xE[i]) = expDivide(expTimes(fl2exp(-U[i, i+1]), (sx[i+1], xE[i+1])), fl2exp(U[i, i]))
    end
    sx, xE;
end

function invBiUexp(U)
    n = size(U)[2]
    (sx, xE) = expBidiagBackSubEn(U)
    (sy, yE) = expInv(expTimes(fl2exp(U[diagind(U)]), (sx, xE)))

    Uinv = zeros(n, n)
    for i in 1:n
        for j in i:n
            Uinvsi, UinvEi = expTimes((sx[i], xE[i]), (sy[j], yE[j]))
            Uinv[i, j] = exp2fl((Uinvsi, UinvEi))
        end
    end
    Uinv;
end


using Test



n = 5000
U = generateTestTriangular(n, 1, Bidiagonal, Float64)

En = one(U)[:, n]
U * (U \ En) ≈ En
x = exp2fl(expBidiagBackSubEn(U))
U * x ≈ En
norm(U * x - En)
U * invBiUexp(U) ≈ I



norm(U * invBiUexp(U) - I)
norm(U * inv(U) - I)


using BenchmarkTools
@btime UinvExp = invBiUexp(U);

@btime Uinv = inv(U);
# UinvEmod = invBidiagUexpmod(U)

