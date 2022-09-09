using BenchmarkTools, LinearAlgebra, BandedMatrices, PyPlot, LaTeXStrings, CSV, DataFrames, Dates


include("backSub.jl")
include("backSubExp.jl")
include("invBiU.jl")
include("invBiUExp.jl")
include("invSymTriSS.jl")
include("invSymTriSSexp.jl")
include("symTriSysSolver.jl")
include("symTriSysSolverexp.jl")
include("genTestMat.jl")



BLAS.get_num_threads()

function time4func(testMat, b, func, solver=false)
    BLAS.set_num_threads(1)
    time = 0.0
    rept = 16
    for i in 1:rept
        if solver
            time += @elapsed func(testMat, b)
        else
            time += @elapsed func(testMat)
        end
    end
    BLAS.set_num_threads(8)
    time / rept;
end





function err4Ainv(A, func)
    Ainv = func(A)
    abs_rightErr = norm(Ainv*A - I)
    abs_leftErr = norm(A*Ainv - I)
    rel_rightErr = abs_rightErr / sqrt(size(A)[1])
    rel_leftErr = abs_leftErr / sqrt(size(A)[1])
    hcat(abs_leftErr, abs_rightErr, rel_leftErr, rel_rightErr);
end


function err4Axb(A, func)
    x = rand(size(A)[1])
    b = A*x
    Ainv = func(A)
    xhat = Ainv * b
    abs_Err = norm(xhat - x)
    rel_Err = abs_Err / norm(x)
    hcat(abs_Err, rel_Err);
end


function err4solver(A, func)
    x = rand(size(A)[1])
    b = A*x
    xhat = func(A, b)
    abs_Err = norm(xhat - x)
    rel_Err = abs_Err / norm(x)
    hcat(abs_Err, rel_Err);
end


function solverCholesky(A, b)
    C = cholesky(Matrix(A))
    U = Bidiagonal(C.U)
    y = U' \ b
    x = U \ y
    x;
end







function CholeskyUtU(U, b)
    y = U' \ b
    x = U \ y
    x;
end


function symTriSS(U, b)
    xhat, yhat = invSymT2SS_U(U, false)
    x = solver4seqTxb(xhat, yhat, b)
    x;
end


function symTriSSexp(U, b)
    (sxh, xhE), (syh, yhE) = invSymT2SSexp_U(U, false)
    x = solver4seqTxbexp((sxh, xhE), (syh, yhE), b)
    x;
end


function seqTxb(func, T, u0, seqlength)
    C = cholesky(Matrix(T))
    U = Bidiagonal(C.U)

    n = length(u0)
    u = zeros(n, seqlength)
    u[:, 1] = u0
    for j in 2:seqlength
        u[:, j] = func(U, u[:, j-1])
    end
    u;
end


function err4seqTxb()
    
end



function getData(funcs, dims, typeM, typeElmt, isTime=true, isLeftRight=true, isTri=false, solver=false)
    if isTime
        data = zeros(length(dims), length(funcs))
    else
        if isLeftRight
            data = zeros(length(dims), 4*length(funcs))
        else
            data = zeros(length(dims), 2*length(funcs))
        end
    end

    for (i, n) in enumerate(dims)
        if isTri
            testMat = generateTestTridiagonal(n, typeM, typeElmt)
            println(typeof(testMat))
        else
            testMat = generateTestTriangular(n, 1, typeM, typeElmt, true)
            println(typeof(testMat))
        end
        b = rand(n)
        for (j, func) in enumerate(funcs)
            if isTime
                BLAS.set_num_threads(1)
                # data[i, j] = time4func(testMat, func, solver)
                data[i, j] = time4func(testMat, b, func, solver)
                # time4func(testMat, func, solver=false, isSeq=false)
                # if solver
                #     data[i, j] = time4func(testMat, func, true, false)
                # else
                #     data[i, j] = time4func(testMat, func, b)
                # end
            else
                BLAS.set_num_threads(8)
                if isLeftRight
                    data[i, 4(j-1)+1:4j] = err4Ainv(testMat, func)
                else
                    if solver
                        data[i, 2(j-1)+1:2j] = err4solver(testMat, func)
                    else
                        data[i, 2(j-1)+1:2j] = err4Axb(testMat, func)
                    end
                end
            end
        end
    end
    BLAS.set_num_threads(8)
    data;
end




function saveData(data, funcs, dims, typeM, typeElmt, isTime=true, isLeftRight=true)
    title_str = String[]
    for func in funcs
        if isTime
            push!(title_str, string(func))
            
        else
            push!(title_str, string(func, "1"))
            push!(title_str, string(func, "2"))
            if isLeftRight
                push!(title_str, string(func, "3"))
                push!(title_str, string(func, "4"))
            end
        end
    end

    df_data = DataFrame(data, :auto)
    rename!(df_data, title_str)
    df_dims = DataFrame(size = dims)
    df = hcat(df_dims, df_data)
    TTIME = string(now())

    if isTime
        CSV.write("/Users/zhiwei_zhou/ic/project/msc_project/formal/invU/data/efficiency/$TTIME$typeM$typeElmt.csv", df)
    else
        CSV.write("/Users/zhiwei_zhou/ic/project/msc_project/formal/invU/data/accuracy/$TTIME$typeM$typeElmt.csv", df)
    end
    TTIME;
end


function saveFig(data, funcs, dims, typeM, typeElmt, TTIME, isTime=true, isLeftRight=true)

    linestylelist = ["--", "-.", ":", "-", ]
    figure()
    grid()
    for (j, func) in enumerate(funcs)
        loglog(dims, data[:, j], linewidth=1.0, linestyle=linestylelist[j], base=2, label=string(func))
    end
    legend()
    if isTime
        xlabel(L"Dimension")
        ylabel("Wall time (s)")
        savefig("/Users/zhiwei_zhou/ic/project/msc_project/formal/invU/figure/efficiency/$TTIME$typeM$typeElmt.png", dpi=150)
        savefig("/Users/zhiwei_zhou/ic/project/msc_project/formal/invU/figure/efficiency/$TTIME$typeM$typeElmt.eps", format="eps")
    else
        xlabel(L"Dimension")
        ylabel("Error")
        savefig("/Users/zhiwei_zhou/ic/project/msc_project/formal/invU/figure/accuracy/$TTIME$typeM$typeElmt.png", dpi=150)
        savefig("/Users/zhiwei_zhou/ic/project/msc_project/formal/invU/figure/accuracy/$TTIME$typeM$typeElmt.eps", format="eps")
    end
    
end




function invTridiag(T)
    qr(T) \ I; 
end

function invBanded(B)
    qr(B) \ I; 
end



function test1()
    isTime=true
    isLeftRight=false

    typeElmt = Float64
    dims = [100*2^i for i = 0:7]
    # dims = [1*2^i for i = 0:3]

    funcs = [inv, invBidiagU, invBiUexp]
    typeM = Bidiagonal
    isTri=false
    data = getData(funcs, dims, typeM, typeElmt, isTime, isLeftRight, isTri)
    TTIME = saveData(data, funcs, dims, typeM, typeElmt, isTime, isLeftRight)
    saveFig(data, funcs, dims, typeM, typeElmt, TTIME, isTime, isLeftRight)

    funcs = [invBanded, invBidiagU, invBiUexp]
    typeM = BandedMatrix
    isTri=false
    data = getData(funcs, dims, typeM, typeElmt, isTime, isLeftRight, isTri)
    TTIME = saveData(data, funcs, dims, typeM, typeElmt, isTime, isLeftRight)
    saveFig(data, funcs, dims, typeM, typeElmt, TTIME, isTime, isLeftRight)


    funcs = [invTridiag, invSymT2SS, invSymT2SSexp]
    typeM = Tridiagonal
    isTri=true
    data = getData(funcs, dims, typeM, typeElmt, isTime, isLeftRight, isTri)
    TTIME = saveData(data, funcs, dims, typeM, typeElmt, isTime, isLeftRight)
    saveFig(data, funcs, dims, typeM, typeElmt, TTIME, isTime, isLeftRight)

    funcs = [invBanded, invSymT2SS, invSymT2SSexp]
    typeM = BandedMatrix
    isTri=true
    data = getData(funcs, dims, typeM, typeElmt, isTime, isLeftRight, isTri)
    TTIME = saveData(data, funcs, dims, typeM, typeElmt, isTime, isLeftRight)
    saveFig(data, funcs, dims, typeM, typeElmt, TTIME, isTime, isLeftRight)

# -----------------------------------------------------------------------------------------------

    isTime=false
    
    typeElmt = Float64
    dims = [100*2^i for i = 0:7]
    # dims = [1*2^i for i = 0:3]

    isLeftRight=false
    funcs = [inv, invBidiagU, invBiUexp]
    typeM = Bidiagonal
    isTri=false
    data = getData(funcs, dims, typeM, typeElmt, isTime, isLeftRight, isTri)
    TTIME = saveData(data, funcs, dims, typeM, typeElmt, isTime, isLeftRight)
    saveFig(data, funcs, dims, typeM, typeElmt, TTIME, isTime, isLeftRight)

    isLeftRight=false
    funcs = [invBanded, invBidiagU, invBiUexp]
    typeM = BandedMatrix
    isTri=false
    data = getData(funcs, dims, typeM, typeElmt, isTime, isLeftRight, isTri)
    TTIME = saveData(data, funcs, dims, typeM, typeElmt, isTime, isLeftRight)
    saveFig(data, funcs, dims, typeM, typeElmt, TTIME, isTime, isLeftRight)



    isLeftRight=true
    funcs = [inv, invBidiagU, invBiUexp]
    typeM = Bidiagonal
    isTri=false
    data = getData(funcs, dims, typeM, typeElmt, isTime, isLeftRight, isTri)
    TTIME = saveData(data, funcs, dims, typeM, typeElmt, isTime, isLeftRight)
    saveFig(data, funcs, dims, typeM, typeElmt, TTIME, isTime, isLeftRight)

    isLeftRight=true
    funcs = [invBanded, invBidiagU, invBiUexp]
    typeM = BandedMatrix
    isTri=false
    data = getData(funcs, dims, typeM, typeElmt, isTime, isLeftRight, isTri)
    TTIME = saveData(data, funcs, dims, typeM, typeElmt, isTime, isLeftRight)
    saveFig(data, funcs, dims, typeM, typeElmt, TTIME, isTime, isLeftRight)



    isLeftRight=false
    funcs = [invTridiag, invSymT2SS, invSymT2SSexp]
    typeM = Tridiagonal
    isTri=true
    data = getData(funcs, dims, typeM, typeElmt, isTime, isLeftRight, isTri)
    TTIME = saveData(data, funcs, dims, typeM, typeElmt, isTime, isLeftRight)
    saveFig(data, funcs, dims, typeM, typeElmt, TTIME, isTime, isLeftRight)

    isLeftRight=false
    funcs = [invBanded, invSymT2SS, invSymT2SSexp]
    typeM = BandedMatrix
    isTri=true
    data = getData(funcs, dims, typeM, typeElmt, isTime, isLeftRight, isTri)
    TTIME = saveData(data, funcs, dims, typeM, typeElmt, isTime, isLeftRight)
    saveFig(data, funcs, dims, typeM, typeElmt, TTIME, isTime, isLeftRight)


    isLeftRight=true
    funcs = [invTridiag, invSymT2SS, invSymT2SSexp]
    typeM = Tridiagonal
    isTri=true
    data = getData(funcs, dims, typeM, typeElmt, isTime, isLeftRight, isTri)
    TTIME = saveData(data, funcs, dims, typeM, typeElmt, isTime, isLeftRight)
    saveFig(data, funcs, dims, typeM, typeElmt, TTIME, isTime, isLeftRight)


    isLeftRight=true
    funcs = [invBanded, invSymT2SS, invSymT2SSexp]
    typeM = BandedMatrix
    isTri=true
    data = getData(funcs, dims, typeM, typeElmt, isTime, isLeftRight, isTri)
    TTIME = saveData(data, funcs, dims, typeM, typeElmt, isTime, isLeftRight)
    saveFig(data, funcs, dims, typeM, typeElmt, TTIME, isTime, isLeftRight)




end

# test1()




# accuracy for [solverCholesky, solverTxbexp]
isTime=false
solver=true
isTri=true
isLeftRight=false

typeElmt = Float64
dims = [100*2^i for i = 0:7]
# dims = [1*2^i for i = 3:5]

funcs = [solverCholesky, solverTxbexp]
typeM = Tridiagonal
isTri=true
data = getData(funcs, dims, typeM, typeElmt, isTime, isLeftRight, isTri, solver)
TTIME = saveData(data, funcs, dims, typeM, typeElmt, isTime, isLeftRight)
# saveFig(data, funcs, dims, typeM, typeElmt, TTIME, isTime, isLeftRight)



[100*2^i for i = 0:7]

for n in [100*2^i for i = 0:7]

    T = Tridiagonal(-ones(n-1), 4ones(n), -ones(n-1))
    u0 = rand(n)
    seqlength = n

    u_UtU = seqTxb(CholeskyUtU, T, u0, seqlength)

    u_SSexp = seqTxb(symTriSSexp, T, u0, seqlength)

    display(norm(u_UtU - u_SSexp))
end

# for random diagonal domiant T
2.519481189322545e-15
5.295133146119094e-15
7.344743001580574e-15
1.1111494812738256e-14
1.70702800397714e-14
2.4321828210194752e-14
3.4970092091237226e-14
5.206321598475564e-14



# for T(-1,4,-1)
2.41261198183003e-14
1.5655785325564103e-13
2.8741597915844873e-13
8.769310262055752e-13
2.8618044076023354e-12
3.693715814742336e-12
1.3379692783891243e-11
6.986714476629926e-11




