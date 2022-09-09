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
    for i in 1:rept
        println((func, i, size(testMat)))
        time += @elapsed func(testMat) # cannot use @belapsed, why?
    end
    time / rept;
end


function timesData4funcs()

    dims = [2^i for i = 2:8]
    dims = [100*2^i for i in 0:10] # for bidiag only
    bw = 1
    typeM = Matrix
    typeM = Bidiagonal
    typeElmt = Float64
    isUpper = true
    rept = 8
    nThreads = 1
    funcs = [inv, invBidiagU, invBiUexp]
    saveData = false
    saveFig = false

    saveData = true
    saveFig = true


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
    xlabel(L"Dimension of $U$"), ylabel("Wall time (s)")
    # ylabel("Wall time (s)")
    grid()
    legend()
    title(L"Time usage for computing $U^{-1}$")
    show()
    if saveFig
        # savefig("./figure/efficiency/$time$typeM$typeElmt.png", dpi=150)
        # savefig("./figure/efficiency/$time$typeM$typeElmt.eps", format="eps")

        savefig("/Users/zhiwei_zhou/ic/project/msc_project/formal/invU/figure/efficiency/$time$typeM$typeElmt.png", dpi=150)
        savefig("/Users/zhiwei_zhou/ic/project/msc_project/formal/invU/figure/efficiency/$time$typeM$typeElmt.eps", format="eps")


    end
    timesData;
end




# timesData4funcs()




n = 100

A = rand(n,n)
T = abs.(Matrix(Tridiagonal(A + A' + n*I)))

U = zero(T)

U[1,1] = sqrt(T[1,1])
U[1,2] = T[1,2]/U[1,1]

for i in 2:n-1
    U[i,i] = sqrt(T[i,i]-U[i-1,i]^2)
    U[i,i+1] = T[i,i+1] / U[i,i]
end
U[n,n] = sqrt(T[n,n] - U[n-1,n]^2)



U

C = cholesky(T)
UU = C.U

U - UU
U' * U - T



