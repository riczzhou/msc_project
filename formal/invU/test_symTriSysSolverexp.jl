using LinearAlgebra, Test

include("invSymTriSSexp.jl")
include("genTestMat.jl")


@testset "Solver of Sym Posi Def Tridiag Linear System (expNum)" begin
    @testset "Float64" begin
        for _ in 1:5
            n = rand(1:1000)
            T = generateTestTridiagonal(n, Matrix, Float64, true)
            b = rand(n)
            x = solveTxbexp(T, b)
            @test T * x ≈ b
        end
    end
end
