using BenchmarkTools, LinearAlgebra, LazyArrays, BandedMatrices, Test
function BackwardSubstitutionV(U, b)
    T = eltype(U)
    n = size(U)[2]
    x = zeros(T, n)
    x[n] = b[n] / U[n, n]
    for i in n-1 : -1 : 1
        sum_ = zero(T)
        for j in i+1 : n
            sum_ += U[i, j] * x[j]
        end
        x[i] = (b[i] - sum_) / U[i, i]
    end
    x
end

function BandedBackwardSubstitutionV(U, b, bw)
    T = eltype(U)
    n = size(U)[2]
    x = zeros(T, n)
    x[n] = b[n] / U[n, n]
    for i in n-1 : -1 : 1
        sum_ = zero(T)
        for j in i+1 : min(n, i+bw)
            sum_ += U[i, j] * x[j]
        end
        x[i] = (b[i] - sum_) / U[i, i]
    end
    x
end

function BandedBackwardSubstitutionM(U, B, bw)
    T = eltype(U)
    n = size(U)[2]
    l = size(B)[2]
    X = zeros(T, n, l)
    for j in l : -1 : 1
        X[:, j] = BandedBackwardSubstitutionV(U, B[:, j], bw)
    end
    X
end

function BackwardSubstitutionM(U, B)
    T = eltype(U)
    n = size(U)[2]
    l = size(B)[2]
    X = zeros(T, n, l)
    for j in l : -1 : 1
        X[:, j] = BackwardSubstitutionV(U, B[:, j])
    end
    X
end

function InverseBidiagonalUpper(U)
    T = eltype(U)
    n = size(U)[2]
    e_n = zeros(T, n)
    e_n[n] = one(T)
    x = BandedBackwardSubstitutionV(U, e_n, 1)
    y = inv.(U[diagind(U)] .* x)
    Uinv = triu(x * y')
    Uinv
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
    D0inv = qr(D0) \ I
    X0 = X[1:r, :]
    Yt0 = invR(D0 * X0)
    Yt[:, 1:r] = Yt0
    Uinv = triu(X * Yt)
    Uinv
end





@testset "Part I" begin
    @testset "Backward substitution" begin
        @testset "For vector b" begin
            for _ in 1:5
                n = rand(1:1000)
                U = big.(triu(rand(n, n))) + I
                b = big.(rand(n))
                x = BackwardSubstitutionV(U, b)
                @test U * x ≈ b
            end
        end
        
        @testset "For banded U and vector b" begin
            for _ in 1:5
                n = rand(1:1000)
                bw = rand(1:n-1)
                U = big.(BandedMatrix(rand(n, n), (0, bw)))
                U += I
                b = big.(rand(n))
                x = BandedBackwardSubstitutionV(U, b, bw)
                @test U * x ≈ b
            end
        end
        
        @testset "For matrix B" begin
            for _ in 1:5
                n = rand(1:1000)
                k = rand(1:10)
                bw = rand(1:n-1)
                U = big.(triu(rand(n, n))) + I
                B = big.(rand(n, k))
                X = BackwardSubstitutionM(U, B)
                @test U * X ≈ B
            end
        end
    
        @testset "For banded U and matrix B" begin
            for _ in 1:5
                n = rand(1:1000)
                k = rand(1:100)
                bw = rand(1:n-1)
                U = big.(BandedMatrix(rand(n, n), (0, bw)))
                U += I
                B = big.(rand(n, k))
                X = BandedBackwardSubstitutionM(U, B, bw)
                @test U * X ≈ B
            end
        end
    end


    @testset "Inverse of Bidiagonal U" begin
        for _ in 1:5
            n = rand(1:1000)
            bw = 1
            U = big.(BandedMatrix(rand(n, n), (0, bw))) + I
            Uinv = InverseBidiagonalUpper(U)
            @test U * Uinv ≈ I
            @test Uinv * U ≈ I
        end
    end

    @testset "Inverse of BlockBidiagonal U" begin
        for _ in 1:5
            n = rand(1:1000)
            bw = rand(1:n-1)
            U = big.(BandedMatrix(rand(n, n), (0, bw))) + I
            Uinv = InverseBlockBiagonalUpper(U, bw)
            @test U * Uinv ≈ I
            @test Uinv * U ≈ I
        end
    end
end

