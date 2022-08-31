using LinearAlgebra, BandedMatrices

function generateTestTriangular(dim, bw, typeM, typeElmt, isUpper=true)
    if bw <= 0 || bw >= dim
        return -1
    end
    
    U = rand(typeElmt, dim, dim) + dim*I
    U = BandedMatrix(U, (0, bw))

    if typeM == Bidiagonal && bw == 1
        U = Bidiagonal(U, :U)
    elseif typeM == UpperTriangular
        U = UpperTriangular(U)
    elseif typeM == Matrix
        U = Matrix(U)
    end
    
    if isUpper
        return U
    else
        return U'
    end
end


function generateTestTridiagonal(dim, typeM, typeElmt, isSymPosiDef=true)
    T = rand(typeElmt, dim, dim) + dim*I
    T = BandedMatrix(T, (1, 1))

    if typeM == Tridiagonal
        T = Tridiagonal(T)
    elseif typeM == Matrix
        T = Matrix(T)
    end

    if isSymPosiDef
        return 0.5(T + T')
    else
        return T
    end
end

