using LinearAlgebra, BandedMatrices, Test


include("backsub.jl")
include("invBiU.jl")
include("genTestMat.jl")

@testset "Inverse of Bidiagonal U" begin
    @testset "Float64" begin
        for _ in 1:5
            n = rand(1:100)
            U = generateTestTriangular(n, 1, Matrix, Float64, true)
            Uinv = invBidiagU(U)
            @test U * Uinv ≈ I
            @test Uinv * U ≈ I
            U = generateTestTriangular(n, 1, Bidiagonal, Float64, true)
            Uinv = invBidiagU(U)
            @test U * Uinv ≈ I
            @test Uinv * U ≈ I
            U = generateTestTriangular(n, 1, BandedMatrix, Float64, true)
            Uinv = invBidiagU(U)
            @test U * Uinv ≈ I
            @test Uinv * U ≈ I
        end
    end

    @testset "BigFloat" begin
        for _ in 1:5
            n = rand(1:1000)
            U = generateTestTriangular(n, 1, Bidiagonal, BigFloat, true)
            Uinv = invBidiagU(U)
            @test U * Uinv ≈ I
            @test Uinv * U ≈ I
        end
    end
end