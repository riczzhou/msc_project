using LinearAlgebra

include("backSub.jl")


function invBidiagU(U)
    T = eltype(U)
    n = size(U)[2]
    e_n = zeros(T, n)
    e_n[n] = one(T)
    x = bandedBackSubVec(U, e_n, 1)
    y = inv.(U[diagind(U)] .* x)
    Uinv = zeros(T, n, n)
    for i in 1:n
        for j in i:n
            Uinv[i, j] = x[i] * y[j]
        end
    end
    # Uinv = triu(x * y')
    Uinv;
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





function invR(A)
    A' * inv(A * A')
end
function invL(A)
    inv(A' * A) * A'
end

function InverseBlockBiagonalUpper(U, bw)
    n = size(U)[2]
    r = rem(n, bw)
    E_k = one(U)[:, n-bw+1:n]
    X = qr(U) \ E_k
    Yt = zero(X')
    for i in n : -bw : r+bw
        Di = U[i-bw+1:i, i-bw+1:i]
        Xi = X[i-bw+1:i, :]
        Yti = qr(Di * Xi) \ I
        Yt[:, i-bw+1:i] = Yti
    end
    D0 = U[1:r, 1:r]
    # D0inv = qr(D0) \ I
    X0 = X[1:r, :]
    Yt0 = invR(D0 * X0)
    Yt[:, 1:r] = Yt0
    Uinv = triu(X * Yt)
    Uinv
end


