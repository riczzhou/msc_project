include("invBiU.jl")


function normCoef(a, b, isReverse=false)
    T = eltype(a)
    n = length(a)
    coef = zeros(T, n)
    if isReverse
        coef[1] = a[1]*b[1]
        for i in 2:n
            coef[i] = coef[i-1] + a[i]*b[i]
        end
    else
        coef[n] = a[n]*b[n]
        for i in n-1:-1:1
            coef[i] = coef[i+1] + a[i]*b[i]
        end
    end
    coef;
end


function generateSymSS(xhat, yhat)
    T = eltype(xhat)
    n = length(xhat)
    S = zeros(T, n, n)
    for i in 1:n
        for j in i:n
            S[i,j] = xhat[i] * yhat[j]
        end
    end
    S + triu(S, 1)';
end


function onepairSymSS(x, y)
    x, x .* normCoef(y,y);
end




function invSymT2SS(T, exactInv=true)
    C = cholesky(Matrix(T))
    U = Bidiagonal(C.U)
    (x, y) = invBidiagU(U,false)
    (xhat, yhat) = onepairSymSS(x, y)
    if exactInv
        S = generateSymSS(xhat, yhat)
        return S;
    end
    xhat, yhat;
end

