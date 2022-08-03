using LinearAlgebra, BandedMatrices, Test


include("backsub.jl")
include("invBiU.jl")


@testset "Backward substitution" begin
    @testset "For vector b" begin
        for _ in 1:5
            n = rand(1:1000)
            U = big.(triu(rand(n, n))) + 100*I
            b = big.(rand(n))
            x = BackSubVec(U, b)
            @test U * x ≈ b
        end
    end

    
    @testset "For banded U and vector b" begin
        for _ in 1:5
            n = rand(1:1000)
            bᵤ = rand(1:n-1)
            U = big.(BandedMatrix(rand(n, n), (0, bᵤ))) + 100*I
            b = big.(rand(n))
            x = BandBackSubVec(U, b, bᵤ)
            @test U * x ≈ b
        end
    end
    
    @testset "For matrix B" begin
        for _ in 1:5
            n = rand(1:1000)
            k = rand(1:10)
            bw = rand(1:n-1)
            U = big.(triu(rand(n, n))) + 100*I
            B = big.(rand(n, k))
            X = BackSubMat(U, B)
            @test U * X ≈ B
        end
    end

    @testset "For banded U and matrix B" begin
        for _ in 1:5
            n = rand(1:1000)
            k = rand(1:100)
            bᵤ = rand(1:n-1)
            U = big.(BandedMatrix(rand(n, n), (0, bᵤ))) + 100*I
            B = big.(rand(n, k))
            X = BandBackSubMat(U, B, bᵤ)
            @test U * X ≈ B
        end
    end
end


@testset "Inverse of Bidiagonal U" begin
    for _ in 1:5
        n = rand(1:1000)
        U = big.(BandedMatrix(rand(n, n), (0, 1))) + 100*I
        Uinv = invBidiagU(U)
        @test U * Uinv ≈ I
        @test Uinv * U ≈ I
    end
end


# @testset "Inverse of BlockBidiagonal U" begin
#     for _ in 1:5
#         n = rand(2:1000)
#         bw = rand(1:n-1)
#         U = big.(BandedMatrix(rand(n, n), (0, bw))) + I
#         Uinv = InverseBlockBiagonalUpper(U, bw)
#         # @time InverseBlockBiagonalUpper(U, bw)
#         # @time qr(U) \ I
#         @test U * Uinv ≈ I
#         @test Uinv * U ≈ I
#     end
# end

