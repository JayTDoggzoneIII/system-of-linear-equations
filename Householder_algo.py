from Matrix_Class import Matrix
from math import fabs
from sqrt_algo import sqrt

EPS = 1e-14

sgn = lambda n: 1 if (n > -EPS) else -1;

smart_mul = lambda A,B,k: [[B[i][j] if (i < k) else sum(A[i][r]*B[r][j] for r in range(k,len(A))) for j in range(len(A))] for i in range(len(A))]

even_more_smart_mul = lambda A,B,k: [[B[i][j] if (i < k or j < k) else sum(A[i][r]*B[r][j] for r in range(k,len(A))) for j in range(len(A))] for i in range(len(A))]

    

def Householder_decomposition(A:Matrix) -> tuple:
    n = len(A)
    R = [[A[i][j] for j in range(n)] for i in range(n)]
    Q = [[0]*n for i in range(n)]   
    for k in range(n-1):
        
        e = [int(i == k) for i in range(k,n)]      
        x = [R[i][k] for i in range(k,n)]
        alpha = -sgn(x[0])*sqrt(sum(x_i*x_i for x_i in x))
        
        u = list(map(lambda p,q: p - alpha * q, x, e))
        norm_u = sqrt(sum(u_i*u_i for u_i in u))
        v = list(map(lambda p: p/norm_u, u))
        
        Q_min = [[int(i == j) - 2*v[i]*v[j] for j in range(n-k)] for i in range(n-k)]
        
        Q_t = [[int(i == j) if (i < k or j < k) else Q_min[i-k][j-k] for j in range(n)] for i in range(n)]
        if (not k):
            Q = Q_t
            R = Matrix(Q_t)*A
        else:
            Q = smart_mul(Q_t,Q,k)
            R = even_more_smart_mul(Q_t,R,k)
        #print(repr(R))
        #print(repr(Q_t))
    return Matrix(Q).T(),Matrix(R);

def solve(A:Matrix, b:list) -> list:
    n = len(b)
    Q,R = Householder_decomposition(A)
    U = Q.T()
    Ub = [sum(U[i][k]*b[k] for k in range(n)) for i in range(n)]
    x = [0]*n
    x[n-1] = Ub[n-1]/R[n-1][n-1]
    for i in reversed(range(n-1)):
        x[i] = (Ub[i] - sum(R[i][j]*x[j] for j in range(i,n)))/R[i][i]
    return x;
