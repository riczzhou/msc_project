function fl2exp(x)
    sx = sign.(x)
    xE = log.(abs.(x))
    sx, xE;
end


function bigfl2exp(x)
    sx = sign.(x)
    xE = log.(abs.(x))
    Float64(sx), Float64(xE);
end


function exp2fl((sx, xE))
    x = sx .* exp.(xE)
    x;
end


function exp2bigfl((sx, xE))
    big.(sx) .* exp.(big.(xE))
end


function expTimes((sx, xE), (sy, yE))
    sxTy = sx .* sy
    xTyE = xE .+ yE
    sxTy, xTyE;
end


function expDivide((sx, xE), (sy, yE))
    sxDy = sx .* sy
    xDyE = xE .- yE
    sxDy, xDyE;
end


function expInv((sx, xE))
    sx, -xE;
end

