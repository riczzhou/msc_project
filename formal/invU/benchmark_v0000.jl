using BenchmarkTools, LinearAlgebra, BandedMatrices, PyPlot, LaTeXStrings, CSV, DataFrames, Dates
include("backsub.jl")
include("invBiU.jl")

include("invBiUExp.jl")


BLAS.get_num_threads()


include("genTestMat.jl")


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
    dims = [2^i for i = 2:8]
    bw = 1
    typeM = Matrix
    typeElmt = Float64
    isUpper = true
    rept = 16
    nThreads = 1
    funcs = [inv, invBidiagU, invBiUexp]
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
        # CSV.write("./data/efficiency/$time$typeM$Float64.csv", df)
        CSV.write("/Users/zhiwei_zhou/ic/project/msc_project/formal/invU/data/efficiency/$time$typeM$Float64.csv", df)
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
        # savefig("./figure/efficiency/$time$typeM$typeElmt.png", dpi=150)
        # savefig("./figure/efficiency/$time$typeM$typeElmt.eps", format="eps")

        savefig("/Users/zhiwei_zhou/ic/project/msc_project/formal/invU/figure/efficiency/$time$typeM$typeElmt.png", dpi=150)
        savefig("/Users/zhiwei_zhou/ic/project/msc_project/formal/invU/figure/efficiency/$time$typeM$typeElmt.eps", format="eps")


    end
    timesData;
end




timesData4funcs()






