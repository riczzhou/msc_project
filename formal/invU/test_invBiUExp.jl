using LinearAlgebra, BandedMatrices, Test


include("backsub.jl")
include("invBiUExp.jl")
include("expNum.jl")
include("genTestMat.jl")


@testset "Inverse of Bidiagonal U (expNum)" begin
    @testset "Back Substitution e_n" begin
        for _ in 1:5
            n = rand(1:1000)
            U = generateTestTriangular(n, 1, Bidiagonal, Float64, true)
            e_n = one(U)[:, n]
            x = exp2fl(expBidiagBackSubEn(U))
            @test U * x ≈ e_n
        end
    end

    @testset "Inverse of BiU" begin
        for _ in 1:5
            n = rand(1:1000)
            U = generateTestTriangular(n, 1, Bidiagonal, Float64, true)
            Uinv = invBiUexp(U)
            @test U * Uinv ≈ I
            @test Uinv * U ≈ I
        end
    end
end


# U * Bidiagonal(invBiUexp(U), :U) ≈ I # banded approx for large n
# norm(U * invBiUexp(U) - I)
# norm(U * Bidiagonal(invBiUexp(U), :U) - I)

