using LinearAlgebra


function BackSubVec(U, b)
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


function BandBackSubVec(U, b, bᵤ)
    T = eltype(U)
    n = size(U)[2]
    x = zeros(T, n)
    x[n] = b[n] / U[n, n]
    for i in n-1 : -1 : 1
        x[i] = b[i]
        for j in i+1 : min(n, i+bᵤ)
            x[i] -= U[i, j] * x[j]
        end
        x[i] /= U[i, i]
    end
    x;
end


function BackSubMat(U, B)
    T = eltype(U)
    n = size(U)[2]
    l = size(B)[2]
    X = zeros(T, n, l)
    for j in l : -1 : 1
        X[:, j] = BackSubVec(U, B[:, j])
    end
    X;
end


function BandBackSubMat(U, B, bᵤ)
    T = eltype(U)
    n = size(U)[2]
    l = size(B)[2]
    X = zeros(T, n, l)
    for j in l : -1 : 1
        X[:, j] = BandBackSubVec(U, B[:, j], bᵤ)
    end
    X;
end

