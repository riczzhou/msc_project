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

