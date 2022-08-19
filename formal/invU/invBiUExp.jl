using LinearAlgebra

include("expNum.jl")


function expBidiagBackSubEn(U)
    n = size(U)[2]
    T = eltype(U)
    sx, xE = zeros(T, n), zeros(T, n)
    (sEn, EnE) = T.((1.0, 0.0))
    (sx[n], xE[n]) = expInv(fl2exp(U[n, n]))
    for i in n-1 : -1 : 1
        (sx[i], xE[i]) = expDivide(expTimes(fl2exp(-U[i, i+1]), (sx[i+1], xE[i+1])), fl2exp(U[i, i]))
    end
    (sx, xE);
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



function getU(n)
    U = rand(n,n) + (n/100)*I
    U -= tril(U, -1) 
    U -= triu(U, 2)
    U
end

n = 1000
U = getU(n)

U = Bidiagonal(rand(n).+100, rand(n-1), :U)


En = one(U)[:, n]
U * (U \ En) ≈ En
x = exp2fl(expBidiagBackSubEn(U))
U * x ≈ En
U * invBiUexp(U) ≈ I


norm(U * invBiUexp(U) - I)
norm(U * inv(U) - I)


using BenchmarkTools
@btime UinvExp = invBiUexp(U)

@btime Uinv = inv(U)
# UinvEmod = invBidiagUexpmod(U)

U = Bidiagonal(rand(n).+100, rand(n-1), :U)