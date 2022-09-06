using LinearAlgebra, Test

include("invSymTriSSexp.jl")
include("genTestMat.jl")



@testset "Inverse of Sym Posi Def Tridiag (expNum)" begin
    @testset "Float64" begin
        for _ in 1:5
            n = rand(1:1000)
            T = generateTestTridiagonal(n, Matrix, Float64, true)
            Tinv = invSymT2SSexp(T)
            @test T * Tinv ≈ I
            @test Tinv * T ≈ I
        end
    end
end