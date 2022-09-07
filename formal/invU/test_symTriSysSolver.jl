using LinearAlgebra, Test

include("symTriSysSolver.jl")
include("genTestMat.jl")


@testset "Solver of Sym Posi Def Tridiag Linear System" begin
    @testset "Float64" begin
        for _ in 1:5
            n = rand(1:50)
            T = generateTestTridiagonal(n, Matrix, Float64, true)
            b = rand(n)
            x = solveTxb(T, b)
            print(norm( T * x - b))
            @test T * x ≈ b
        end
    end

    @testset "BigFloat" begin
        for _ in 1:5
            n = rand(1:100)
            T = generateTestTridiagonal(n, Matrix, BigFloat, true)
            b = rand(n)
            x = solveTxb(T, b)
            @test T * x ≈ b
        end
    end
end

