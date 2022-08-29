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

# @testset "Inverse of Bidiagonal U" begin
#     for _ in 1:1
#         n = rand(1:1000)
#         bw = 1
#         U = big.(BandedMatrix(rand(n, n), (0, bw))) + I
#         Uinv = InverseBidiagonalUpper(U)
#         @test U * Uinv ≈ I
#         @test Uinv * U ≈ I
#     end
# end