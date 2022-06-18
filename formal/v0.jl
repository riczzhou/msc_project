
function BackwardSubstitution(U, b)
    T = eltype(U)
    n = size(U)[2]
    x = zeros(T, n)
    x[n] = b[n] / U[n, n]
    for i in n-1 : -1 : 1
        sum_ = zero(T)
        for j in i+1 : n
            sum_ += U[i, j] * x[j]
        end
        x[i] = (b[i] - sum_) / U[i, i]
    end
    x
end


function BandedBackwardSubstitution(U, b, bw)
    T = eltype(U)
    n = size(U)[2]
    x = zeros(T, n)
    x[n] = b[n] / U[n, n]
    for i in n-1 : -1 : 1
        sum_ = zero(T)
        for j in i+1 : min(n, i+bw)
            sum_ += U[i, j] * x[j]
        end
        x[i] = (b[i] - sum_) / U[i, i]
    end
    x
end


function TEST_BackwardSubstitution()
    for _ in 1:10
        n = rand(1:2000)
        U = big.(triu(rand(n, n)))
        U += I
        b = big.(rand(n))
        x = BackwardSubstitution(U, b)
        println(U * x ≈ b)
    end
end

function TEST_BandedBackwardSubstitution()
    for _ in 1:10
        n = rand(1:2000)
        bw = rand(1:n)
        U = big.(BandedMatrix(rand(n, n), (0, bw)))
        U += I
        b = big.(rand(n))
        x = BandedBackwardSubstitution(U, b, bw)
        println(U * x ≈ b)
    end
end


TEST_BackwardSubstitution()
TEST_BandedBackwardSubstitution()



function InvseBidiagonalUpper(U)
    T = eltype(U)
    n = size(U)[2]
    e_n = zeros(T, n)
    e_n[n] = one(T)
    x = BandedBackwardSubstitution(U, e_n, 1)
    y = inv.(U[diagind(U)] .* x)
    Uinv = triu(x * y')
    Uinv
end



function TEST_InvseBidiagonalUpper()
    for _ in 1:10
        n = rand(1:2000)
        # bw = rand(1:n)
        bw = 1
        U = big.(BandedMatrix(rand(n, n), (0, bw)))
        U += I
        Uinv = InvseBidiagonalUpper(U)
        println(U * Uinv ≈ I)
    end
end

TEST_InvseBidiagonalUpper()





