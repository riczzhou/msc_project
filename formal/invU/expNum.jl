function fl2exp(x)
    sx = sign.(x)
    xE = log.(abs.(x))
    (sx, xE);
end


function exp2fl((sx, xE))
    x = sx .* exp.(xE);
end


function expTimes((sx, xE), (sy, yE))
    sxTy = sx .* sy
    xTyE = xE .+ yE
    (sxTy, xTyE);
end


function expDivide((sx, xE), (sy, yE))
    sxDy = sx .* sy
    xDyE = xE .- yE
    (sxDy, xDyE);
end


function expInv((sx, xE))
    (sx, -xE);
end


n = 100
x = rand(n)
y = rand(n)

exp2fl(fl2exp(x)) ≈ x
exp2fl(expTimes(fl2exp(x), fl2exp(y))) ≈ x .* y
exp2fl(expDivide(fl2exp(x), fl2exp(y))) ≈ x ./ y
exp2fl(expInv(fl2exp(x))) ≈ 1 ./ x
