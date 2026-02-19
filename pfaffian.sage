R = PolynomialRing(QQ, 28, 'x')
vars = R.gens()

X = matrix(R, 8, 8)

k = 0
for i in range(8):
    for j in range(i):
        X[i,j] = vars[k]
        X[j,i] = -vars[k]
        k += 1


J1 = matrix(R, [[0,-1],[1,0]])
J = block_diagonal_matrix(J1, J1, J1, J1)


def jordan(x,y):
    return (x*J*y + y*J*x)/2

def trd(x):
    return (1/2)*(x*J).trace()

def pairing(x,y):
    return trd(jordan(x,y))

def quartic(X):
    xx = jordan(X,X)
    return (1/24)*(
        trd(X)^4
        - 6*trd(X)^2*pairing(X,X)
        + 3*pairing(X,X)^2
        + 8*trd(X)*pairing(xx,X)
        - 6*pairing(xx,xx)
    )

Q = quartic(X)
P = X.pfaffian()

print((Q - P).is_zero())
