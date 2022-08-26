using BenchmarkTools, LinearAlgebra, BandedMatrices, PyPlot, LaTeXStrings, CSV, DataFrames, Dates
include("backsub.jl")
include("invBiU.jl")

BLAS.get_num_threads()

function generateTestTriangular(dim, bw, typeM, typeElmt, isUpper=true)
    if bw <= 0
        return -1
    end
    
    U = rand(typeElmt, dim, dim) + dim*I

    if typeM == Matrix
        U = triu(U) - triu(U, bw+1)
    elseif typeM == Bidiagonal
        bw == 1 ? U = Bidiagonal(U, :U) : return -1
    elseif typeM == BandedMatrix
        bw <= n-1 ? U = BandedMatrix(U, (0, bw)) : return -1
    elseif typeM == UpperTriangular
        bw == n-1 ? U = UpperTriangular(U) : return -1
    else
        return -1
    end

    if isUpper
        return U
    else
        return U'
    end
end


function timeAvg4func(func, testMat, rept, nThreads)
    println(func)
    BLAS.set_num_threads(nThreads)
    time = 0.0
    for _ in 1:rept
        time += @elapsed func(testMat) # cannot use @belapsed, why?
    end
    time / rept;
end


function timesData4funcs()
    # M = generateTestTriangular(n, bw, typeM, typeElmt, isUpper)
    # M = generateTestTridiagonal(n, typeM, typeElmt, isSymPosiDef)
    dims = [2^i for i = 4:9]
    bw = 1
    typeM = Matrix
    typeElmt = Float64
    isUpper = true
    rept = 10
    nThreads = 1
    funcs = [inv, invBidiagU]
    saveData = false
    saveFig = false


    timesData = zeros(length(dims), length(funcs))

    BLAS.set_num_threads(nThreads)
    for (j, func) in enumerate(funcs)
        for (i, n) in enumerate(dims)
            testMat = generateTestTriangular(n, bw, typeM, typeElmt, isUpper)
            timesData[i, j] = timeAvg4func(func, testMat, rept, nThreads)
        end
    end
    BLAS.set_num_threads(8)

    if saveData
        df_data = DataFrame(timesData, :auto)
        rename!(df_data, string.(funcs))
        df_dims = DataFrame(size = dims)
        df = hcat(df_dims, df_data)
        time = string(now())
        CSV.write("./data/efficiency/$time$typeM$Float64.csv", df)
        # CSV.write("./data/accuracy/$time$typeM$Float64.csv", df)
    end

    linestylelist = ["--", "-.", ":", "-", ]
    figure()
    for (j, func) in enumerate(funcs)
        # loglog(dims, timesData[:, j], color=colorlist[j], linewidth=1.0, linestyle=linestylelist[j], base=2, label=string(func))
        loglog(dims, timesData[:, j], linewidth=1.0, linestyle=linestylelist[j], base=2, label=string(func))
    end
    xlabel(L"Dimension of $U$, $N$"), ylabel("Wall time (s)")
    # ylabel("Wall time (s)")
    grid()
    legend()
    title(L"Time usage for computing $U^{-1} \in \mathbb{R}^{N \times N}$")
    show()
    if saveFig
        savefig("./figure/efficiency/$time$typeM$typeElmt.png", dpi=150)
        savefig("./figure/efficiency/$time$typeM$typeElmt.eps", format="eps")
    end
    timesData;
end
timesData4funcs()






