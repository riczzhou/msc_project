function backSubVec(U, b)
    T, n = eltype(U), size(U)[2]
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


function invBidiagU(U, exactInv=true)
    T = eltype(U)
    n = size(U)[2]
    e_n = one(U)[:,n]
    x = bandedBackSubVec(U, e_n, 1)
    y = inv.(U[diagind(U)] .* x)
    if exactInv
        Uinv = zeros(T, n, n)
        for i in 1:n
            for j in i:n
                Uinv[i, j] = x[i] * y[j]
            end
        end
        return Uinv;
    end
    x, y;
end