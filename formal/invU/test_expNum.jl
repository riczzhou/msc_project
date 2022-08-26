using Test

include("expNum.jl")

@testset "Exponential Number" begin
    @testset "Transformation" begin
        for _ in 1:5
            n = rand(1:10000)
            x = rand(n)
            @test exp2fl(fl2exp(x)) ≈ x
        end
    end

    @testset "Multiplication" begin
        for _ in 1:5
            n = rand(1:10000)
            x = rand(n)
            y = rand(n)
            @test exp2fl(expTimes(fl2exp(x), fl2exp(y))) ≈ x .* y
        end
    end

    @testset "Division" begin
        for _ in 1:5
            n = rand(1:10000)
            x = rand(n)
            y = rand(n)
            @test exp2fl(expDivide(fl2exp(x), fl2exp(y))) ≈ x ./ y
        end
    end

    @testset "Inverse" begin
        for _ in 1:5
            n = rand(1:10000)
            x = rand(n)
            y = rand(n)
            @test exp2fl(expInv(fl2exp(x))) ≈ 1.0 ./ x
        end
    end
end