using LinearAlgebra

include("backSub.jl")


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
    # D0inv = qr(D0) \ I
    X0 = X[1:r, :]
    Yt0 = invR(D0 * X0)
    Yt[:, 1:r] = Yt0
    Uinv = triu(X * Yt)
    Uinv
end


include("genTestMat.jl")


n = 100
bw = 3
U = generateTestTriangular(n, bw, BandedMatrix, BigFloat)
# U = generateTestTriangular(n, bw, Bidiagonal, Float64)

# U * invBidiagU(U) ≈ I
U * InverseBlockBiagonalUpper(U, bw) ≈ I

# invBidiagU(U)


# @testset "Inverse of BlockBidiagonal U" begin
#     for _ in 1:1
#         n = 1000
#         bw = 97
#         U = big.(BandedMatrix(rand(n, n), (0, bw))) + I
#         Uinv = InverseBlockBiagonalUpper(U, bw)
#         # @time InverseBlockBiagonalUpper(U, bw)
#         # @time qr(U) \ I
#         @test U * Uinv ≈ I
#         @test Uinv * U ≈ I
#     end
# end