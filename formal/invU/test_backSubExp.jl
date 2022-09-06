using LinearAlgebra, BandedMatrices, Test

include("backSubExp.jl")
include("genTestMat.jl")

n = 100
bw = 32
U = generateTestTriangular(n, n-1, UpperTriangular, Float64)
U * exp2fl(backSubMatExp(U, one(U))) ≈ I


U = generateTestTriangular(n, bw, BandedMatrix, Float64)
U * exp2fl(bandedBackSubMatExp(U, one(U), bw)) ≈ I







bw = 32
U = generateTestTriangular(n, bw, BandedMatrix, Float64)

b = rand(n)
xx = exp2fl(backSubVecExp(U, b))
norm(U * xx - b)

yy = exp2fl(bandedBackSubVecExp(U, b, bw))
norm(U * yy - b)


@testset "Backward Substitution (expNum)" begin
    @testset "U \\ b (expNum)" begin
        for _ = 1:5
            n = rand(1:1000)
            bw = rand(1:n-1)
            U = generateTestTriangular(n, bw, Matrix, Float64, true)
            b = rand(n)
            x = exp2fl(backSubVecExp(U, b))
            @test U * x ≈ b
        end
    end

    @testset "banded U \\ b (expNum)" begin
        for _ = 1:5
            n = rand(1:1000)
            bw = rand(1:n-1)
            U = generateTestTriangular(n, bw, BandedMatrix, Float64, true)
            b = rand(n)
            x = exp2fl(bandedBackSubVecExp(U, b, bw))
            @test U * x ≈ b
        end
    end

    @testset "U \\ B (expNum)" begin
        for _ = 1:5
            n = rand(1:1000)
            k = rand(1:10)
            bw = rand(1:n-1)
            U = generateTestTriangular(n, bw, Matrix, Float64, true)
            B = rand(n, k)
            X = exp2fl(backSubMatExp(U, B))
            @test U * X ≈ B
        end
    end

    @testset "banded U \\ B (expNum)" begin
        for _ = 1:5
            n = rand(1:1000)
            k = rand(1:10)
            bw = rand(1:n-1)
            U = generateTestTriangular(n, bw, BandedMatrix, Float64, true)
            B = rand(n, k)
            X = exp2fl(bandedBackSubMatExp(U, B, bw))
            @test U * X ≈ B
        end
    end
end

