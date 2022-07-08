using BenchmarkTools, LinearAlgebra, LazyArrays, BandedMatrices


function backSub(U, b)
    T = eltype(U)
    n = size(U)[2]
    x = zeros(T, n)
    x[end] = b[end] / U[end, end]
    for i in n-1 : -1 : 1
        sum_ = 0.0
        for j in i+1 : n
            sum_ += U[i, j] * x[j]
        end
        x[i] = (b[i] - sum_)/U[i, i]
    end
    x
end

function backSubBanded(U, b, bw=1)
    T = eltype(U)
    n = size(U)[2]
    x = zeros(T, n)
    x[end] = b[end] / U[end, end]
    for i in n-1 : -1 : 1
        sum_ = zero(T)
        for j in i+1 : min(n, i+bw)
            sum_ += U[i, j] * x[j]
        end
        x[i] = (b[i] - sum_) / U[i, i]

    end
    x
end


function invBiU(U)
    T = eltype(U)
    n = size(U)[2]
    e_n = zeros(T, n)
    e_n[n] = one(T)
    # x = U \ e_n
    x = backSubBanded(U, e_n, 1)
    y = 1 ./ x
    Uinv = triu(x * y')
    Uinv
end



function invU_1_0(U)
    T = eltype(U)
    n = size(U)[2]
    e_n = zeros(T, n)
    e_n[n] = one(T)
    # x = U \ e_n
    x = backSubBanded(U, b, 1)
    y = 1 ./ (U[diagind(U)] .* x)
    # y = 1 ./ x
    Uinv = triu(x * y')
    Uinv
end



n = 100
bw = 1


U = big.(Bidiagonal( 1 .+ rand(n), rand(n-1), :U))

b = big.(rand(n))

norm((U \ b) - backSub(U, b))

norm((U \ b) - backSubBanded(U, b, 1))

typeof(1 ./ backSubBanded(U, b, 1))


norm(invU_1_0(U) - inv(U))




A = big.(rand(n, n))
A = triu(A) - triu(A, 2)
# A += I

A[diagind(A)] .= big.(ones(n))
# A += I
A
b = big.(zeros(n))
b[n] = big.(one(1))
b
x = A \ b

x = backSubBanded(A, b, 1)

A * x ≈ b

A[diagind(A)]
y = inv.(A[diagind(A)] .* x)
y = 1 ./ x
norm(triu(x * y') - inv(U))




triu(x * y')
inv(U)