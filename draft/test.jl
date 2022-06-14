using LinearAlgebra, BandedMatrices
n = 4
bw = 2
U = BandedMatrix(rand(n, n), (0, bw))
for i in 1:n
    U[i, i] = 1.0
end
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
Y1 = -inv(D1*tril(X1, -1) + B1*X2)

Yt = [Y1 Y2]
triu(X * Yt) - inv(U)
norm(triu(X * Yt) - inv(U))

