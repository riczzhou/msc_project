using LinearAlgebra

include("backSub.jl")


function invBidiagU(U, exactInv=true)
    T = eltype(U)
    n = size(U)[2]
    e_n = zeros(T, n)
    e_n[n] = one(T)
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


function invBidiagUxy(U)
    T = eltype(U)
    n = size(U)[2]
    e_n = zeros(T, n)
    e_n[n] = one(T)
    x = BandBackSubVec(U, e_n, 1)
    y = inv.(U[diagind(U)] .* x)
    Uinv = zeros(T, n, n)
    for i in 1:n
        for j in i:n
            Uinv[i, j] = x[i] * y[j]
        end
    end
    # Uinv = triu(x * y')
    (Uinv, x, y);
end

