function normCoefexp((sa, aE), (sb, bE), isReverse=false)
    n = length(sa)
    scoef, coefE = zeros(n), zeros(n)
    if isReverse
        scoef[1], coefE[1] = expTimes((sa[1], aE[1]), (sb[1], bE[1]))
        for i in 2:n
            fl1 = exp2bigfl((scoef[i-1], coefE[i-1]))
            fl2 = exp2bigfl(expTimes((sa[i], aE[i]), (sb[i], bE[i])))
            scoef[i], coefE[i] = bigfl2exp(fl1+fl2)
        end
    else
        scoef[n], coefE[n] = expTimes((sa[n], aE[n]), (sb[n], bE[n]))
        for i in n-1:-1:1
            fl1 = exp2bigfl((scoef[i+1], coefE[i+1]))
            fl2 = exp2bigfl(expTimes((sa[i], aE[i]), (sb[i], bE[i])))
            scoef[i], coefE[i] = bigfl2exp(fl1+fl2)
        end
    end
    scoef, coefE;
end

function onepairSymSSexp((sx, xE), (sy, yE))
    (sx, xE), expTimes((sx, xE), normCoefexp((sy, yE), (sy, yE)));
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