using LinearAlgebra, BandedMatrices, Test

include("backSub.jl")
include("genTestMat.jl")


@testset "Backward Substitution" begin
    @testset "U \\ b" begin
        for _ = 1:5
            n = rand(1:1000)
            bw = rand(1:n-1)
            U = generateTestTriangular(n, bw, Matrix, Float64, true)
            b = rand(n)
            x = backSubVec(U, b)
            @test U * x ≈ b
            U = generateTestTriangular(n, bw, Matrix, BigFloat, true)
            b = big.(rand(n))
            x = backSubVec(U, b)
            @test U * x ≈ b
        end
    end

    @testset "banded U \\ b" begin
        for _ = 1:5
            n = rand(1:1000)
            bw = rand(1:n-1)
            U = generateTestTriangular(n, bw, BandedMatrix, Float64, true)
            b = rand(n)
            x = bandedBackSubVec(U, b, bw)
            @test U * x ≈ b
            U = generateTestTriangular(n, bw, BandedMatrix, BigFloat, true)
            b = big.(rand(n))
            x = bandedBackSubVec(U, b, bw)
            @test U * x ≈ b
        end
    end
    
    @testset "U \\ B" begin
        for _ = 1:5
            n = rand(1:1000)
            k = rand(1:10)
            bw = rand(1:n-1)
            U = generateTestTriangular(n, bw, Matrix, Float64, true)
            B = rand(n, k)
            X = backSubMat(U, B)
            @test U * X ≈ B
            U = generateTestTriangular(n, bw, Matrix, BigFloat, true)
            B = big.(rand(n, k))
            X = backSubMat(U, B)
            @test U * X ≈ B
        end
    end

    @testset "banded U \\ B" begin
        for _ = 1:5
            n = rand(1:1000)
            k = rand(1:10)
            bw = rand(1:n-1)
            U = generateTestTriangular(n, bw, BandedMatrix, Float64, true)
            B = rand(n, k)
            X = backSubMat(U, B)
            @test U * X ≈ B
            U = generateTestTriangular(n, bw, BandedMatrix, BigFloat, true)
            B = big.(rand(n, k))
            X = backSubMat(U, B)
            @test U * X ≈ B
        end
    end
end

