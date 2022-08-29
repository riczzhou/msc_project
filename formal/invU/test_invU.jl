using LinearAlgebra, BandedMatrices, Test


include("backsub.jl")
include("invBiU.jl")



@testset "Inverse of Bidiagonal U" begin
    for _ in 1:5
        n = rand(1:1000)
        U = big.(BandedMatrix(rand(n, n), (0, 1))) + 100*I
        Uinv = invBidiagU(U)
        @test U * Uinv ≈ I
        @test Uinv * U ≈ I
    end
end



