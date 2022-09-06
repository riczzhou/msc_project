using LinearAlgebra, Test

include("invSymTriSS.jl")
include("genTestMat.jl")



@testset "Inverse of Sym Posi Def Tridiag" begin
    @testset "Float64" begin
        for _ in 1:5
            n = rand(1:100)
            T = generateTestTridiagonal(n, Matrix, Float64, true)
            Tinv = invSymT2SS(T)
            @test T * Tinv ≈ I
            @test Tinv * T ≈ I
        end
    end

    @testset "BigFloat" begin
        for _ in 1:5
            n = rand(1:1000)
            T = generateTestTridiagonal(n, Matrix, BigFloat, true)
            Tinv = invSymT2SS(T)
            @test T * Tinv ≈ I
            @test Tinv * T ≈ I
        end
    end
end