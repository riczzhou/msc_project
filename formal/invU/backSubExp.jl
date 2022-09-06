using LinearAlgebra

include("expNum.jl")
include("genTestMat.jl")


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


function backSubMatExp(U, B)
    T = eltype(U)
    n = size(U)[2]
    l = size(B)[2]
    sX, XE = zeros(T, n, l), zeros(T, n, l)

    for j in l : -1 : 1
        (sX[:, j], XE[:, j]) = backSubVecExp(U, B[:, j])
    end
    sX, XE;
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


function bandedBackSubMatExp(U, B, bw)
    T = eltype(U)
    n = size(U)[2]
    l = size(B)[2]
    sX, XE = zeros(T, n, l), zeros(T, n, l)
    for j in l : -1 : 1
        (sX[:, j], XE[:, j]) = bandedBackSubVecExp(U, B[:, j], bw)
    end
    sX, XE;
end
