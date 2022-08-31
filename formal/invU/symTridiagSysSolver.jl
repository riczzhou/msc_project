include("invBiUexp.jl")
include("genTestMat.jl")

dim = 100
typeM = Tridiagonal
typeElmt = Float64
isSymPosiDef=true

T = generateTestTridiagonal(dim, typeM, typeElmt, isSymPosiDef)
invBiUexp(T)




# function normCoef(u, v, isReverse=True)
#     T = eltype(u)
#     n = length(u)
#     Coef = zeros(T, n)
#     Coef[1] = u[n]*v[n]
#     for i in n-1:-1:1
#         coef[i] = coef[i+1] + u[i]*v[i]
#     end
#     coef;
# end



# a = [1,2,3]
# reverse(a)