using LinearAlgebra, BandedMatrices

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


function generateTestTridiagonal(dim, typeM, typeElmt, isSymPosiDef=true)
    T = rand(typeElmt, dim, dim) + dim*I

    if typeM == Matrix
        T = triu(T, -1) - triu(T, 2)
    elseif typeM == Tridiagonal
        T = Tridiagonal(T)
    elseif typeM == BandedMatrix
        T = BandedMatrix(T, (1, 1))
    else
        return -1
    end

    if isSymPosiDef
        return 0.5(T + T')
    else
        return T
    end
end

