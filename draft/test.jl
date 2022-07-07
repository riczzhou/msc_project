using LinearAlgebra, BandedMatrices
n = 4
bw = 2
U = BandedMatrix(rand(n, n), (0, bw))
# for i in 1:n
#     U[i, i] = 1.0
# end
display(U)

D1 = U[1:2, 1:2]
D2 = U[3:4, 3:4]
B1 = U[1:2, 3:4]

E2 = [zeros(2, 2) ; I]
E1 = [I ; zeros(2, 2)]

X = U \ E2
X1 = X[1:2, :]
X2 = X[3:4, :]

Y2 = [1 0 ; 0 1]
Y1 = inv(X1) * inv(D1)

Yt = [Y1 Y2]
triu(X * Yt) - inv(U)
norm(triu(X * Yt) - inv(U))


n = 2*2
bw = 2
U = BandedMatrix(rand(n, n), (0, bw))
# for i in 1:n
#     U[i, i] = 1.0
# end
display(U)


E2 = zeros(n, 2)

E2[end-1:end, :] = [1 0 ; 0 1]
X = U \ E2
Y = zeros(n, 2)
Y[n-1:n, :] = [1 0 ; 0 1]


a = 0
println(X)
for i in n-2 : -2 : 1
    Xi = X[i: i+1, :]
    print(Xi)
end 
a


# A = similar(X)

[i for i in n-2 : -2 : 0]

using LinearAlgebra
Threads.nthreads()

BLAS.get_num_threads()

BLAS.set_num_threads(8)
BLAS.get_num_threads()



include("testforimport.jl")

testforimport(1,2)

