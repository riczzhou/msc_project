using LinearAlgebra

function backSubVec(U, b)
    T = eltype(U)
    n = size(U)[2]
    x = zeros(T, n)
    x[n] = b[n] / U[n, n]
    for i in n-1 : -1 : 1
        x[i] = b[i]
        for j in i+1 : n
            x[i] -= U[i, j] * x[j]
        end
        x[i] /= U[i, i]
    end
    x;
end


function bandedBackSubVec(U, b, bw)
    T = eltype(U)
    n = size(U)[2]
    x = zeros(T, n)
    x[n] = b[n] / U[n, n]
    for i in n-1 : -1 : 1
        x[i] = b[i]
        for j in i+1 : min(n, i+bw)
            x[i] -= U[i, j] * x[j]
        end
        x[i] /= U[i, i]
    end
    x;
end


function backSubMat(U, B)
    T = eltype(U)
    n = size(U)[2]
    l = size(B)[2]
    X = zeros(T, n, l)
    for j in l : -1 : 1
        X[:, j] = backSubVec(U, B[:, j])
    end
    X;
end


function bandedBackSubMat(U, B, bw)
    T = eltype(U)
    n = size(U)[2]
    l = size(B)[2]
    X = zeros(T, n, l)
    for j in l : -1 : 1
        X[:, j] = bandedBackSubVec(U, B[:, j], bw)
    end
    X;
end

