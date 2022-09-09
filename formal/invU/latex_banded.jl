function invR(A)
    A' * inv(A * A')
end
function invL(A)
    inv(A' * A) * A'
end

function InverseBlockBiagonalUpper(U, bw) # have problem need note
    n = size(U)[2]
    r = rem(n, bw)
    T = eltype(U)
    E_k = one(T, U)[:, n-bw+1:n]
    X = qr(U) \ E_k
    Yt = zero(T, X')
    for i in n : -bw : r+bw
        Di = U[i-bw+1:i, i-bw+1:i]
        Xi = X[i-bw+1:i, :]
        Yti = qr(Di * Xi) \ I
        Yt[:, i-bw+1:i] = Yti
    end
    D0 = U[1:r, 1:r]
    X0 = X[1:r, :]
    Yt0 = invR(D0 * X0)
    Yt[:, 1:r] = Yt0
    Uinv = triu(X * Yt)
    Uinv
end