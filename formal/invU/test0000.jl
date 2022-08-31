using LinearAlgebra

include("backSub.jl")




function invBidiagUxy(U)
    T = eltype(U)
    n = size(U)[2]
    e_n = zeros(T, n)
    e_n[n] = one(T)
    x = bandedBackSubVec(U, e_n, 1)
    y = inv.(U[diagind(U)] .* x)
    (x, y);
end


function innerCoef(u, v)
    T = eltype(u)
    n = length(u)
    coef = zeros(T, n)
    coef[n] = u[n]*v[n]
    for i in n-1:-1:1
        coef[i] = coef[i+1] + u[i]*v[i]
    end
    coef;
end


function symSemisep(x, y)
    x̃ = x
    ỹ = x .* innerCoef(y,y)
    (x̃, ỹ);
end


function generateSymSS(x, y)
    # tril(ỹ * x̃') + triu(x̃ * ỹ', 1)
    T = eltype(x)
    n = length(x)
    SS = zeros(T, n, n)
    for i in 1:n
        for j in 1:n
            if i <= j
                SS[i,j] = y[j]*x[i]
            else
                SS[i,j] = y[i]*x[j]
            end
        end
    end
    SS
end


function invSymTridiagSS(T)
    C = cholesky(T)
    U = C.U
    (x, y) = invBidiagUxy(U)
    (x̃, ỹ) = symSemisep(x, y)
    # SS = tril(ỹ * x̃') + triu(x̃ * ỹ', 1)
    SS = generateSymSS(x̃, ỹ)
    # println(norm(SS-inv(T)))
    SS
end




n = 1000
A = big.(Tridiagonal(rand(n,n)))
T = Matrix(A + A' + 100*I)
# @btime invSymTridiagSS(T);
# @btime inv(T);
# invSymTridiagSS(T)




Tinv = invSymTridiagSS(T)

Tinv * T ≈ I


using BenchmarkTools



