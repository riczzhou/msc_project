include("invBiUexp.jl")
include("genTestMat.jl")

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








function exp2bigfl((sx, xE))
    big.(sx) .* exp.(big.(xE))
end


function bigfl2exp(x)
    sx = sign.(x)
    xE = log.(abs.(x))
    Float64(sx), Float64(xE);
end


function normCoefexp((sa, aE), (sb, bE), isReverse=false)
    n = length(sa)
    scoef, coefE = zeros(n), zeros(n)
    if isReverse
        scoef[1], coefE[1] = expTimes((sa[1], aE[1]), (sb[1], bE[1]))
        for i in 2:n
            # fl1 = exp2fl((scoef[i-1], coefE[i-1]))
            # fl2 = exp2fl(expTimes((sa[i], aE[i]), (sb[i], bE[i])))
            # scoef[i], coefE[i] = fl2exp(fl1+fl2)

            fl1 = exp2bigfl((scoef[i-1], coefE[i-1]))
            fl2 = exp2bigfl(expTimes((sa[i], aE[i]), (sb[i], bE[i])))
            scoef[i], coefE[i] = bigfl2exp(fl1+fl2)
        end
    else
        scoef[n], coefE[n] = expTimes((sa[n], aE[n]), (sb[n], bE[n]))
        for i in n-1:-1:1
            # fl1 = exp2fl((scoef[i+1], coefE[i+1]))
            # fl2 = exp2fl(expTimes((sa[i], aE[i]), (sb[i], bE[i])))
            # scoef[i], coefE[i] = fl2exp(fl1+fl2)

            fl1 = exp2bigfl((scoef[i+1], coefE[i+1]))
            fl2 = exp2bigfl(expTimes((sa[i], aE[i]), (sb[i], bE[i])))
            scoef[i], coefE[i] = bigfl2exp(fl1+fl2)
        end
    end
    # println(coefE)
    # println()

    scoef, coefE;
end





a = [1,2,3]
normCoef(a, a, true)
scoef, coefE = normCoefexp(fl2exp(a), fl2exp(a), true)
exp2fl((scoef, coefE))

normCoef(a, a)
scoef, coefE = normCoefexp(fl2exp(a), fl2exp(a))
exp2fl((scoef, coefE))




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



function generateSymSSexp((sxh, xhE), (syh, yhE))
    n = length(sxh)
    S = zeros(n, n)
    for i in 1:n
        for j in i:n
            sSi, SiE = expTimes((sxh[i], xhE[i]), (syh[j], yhE[j]))
            S[i, j] = exp2fl((sSi, SiE))
        end
    end
    S + triu(S, 1)';
end



generateSymSS(a,a)
generateSymSSexp(fl2exp(a),fl2exp(a))



function onepairSymSS(x, y)
    x, x .* normCoef(y,y);
end



function onepairSymSSexp((sx, xE), (sy, yE))
    (sx, xE), expTimes((sx, xE), normCoefexp((sy, yE), (sy, yE)));
end



onepairSymSS(a,a)
(sx, xE), (sy, yE) = onepairSymSSexp(fl2exp(a),fl2exp(a))
exp2fl((sx, xE)), exp2fl((sy, yE))


function invSymT2SS(T)
    C = cholesky(T)
    U = C.U
    U = Bidiagonal(C.U)
    (x, y) = invBidiagU(U,false)
    (xhat, yhat) = onepairSymSS(x, y)
    SS = generateSymSS(xhat, yhat)
    SS
end



function invSymT2SSexp(T)
    C = cholesky(T)
    U = C.U
    U = Bidiagonal(C.U)
    (sx, xE), (sy, yE) = invBiUexp(U,false)
    # println((sx, xE), (sy, yE))
    (sxh, xhE), (syh, yhE) = onepairSymSSexp((sx, xE), (sy, yE))
    # println()
    # println(xhE)
    # println()
    # println(yhE)
    SS = generateSymSSexp((sxh, xhE), (syh, yhE))
    SS
end






n = 50000
T = Tridiagonal(-ones(n-1), 4ones(n), -ones(n-1));
T = Matrix(T);
# T = big.(Matrix(T));
# Tinv = invSymT2SS(T)
# norm(T * Tinv - I)


Tinv = invSymT2SSexp(T);
norm(T * Tinv - I)


# exp(Inf)



# dim = 100
# typeM = Tridiagonal
# typeElmt = Float64
# isSymPosiDef=true

# T = generateTestTridiagonal(dim, typeM, typeElmt, isSymPosiDef)
# invBiUexp(T)


# invBidiagU(U, false)



Float64(big(big(10)^300))
big(big(10)^30)