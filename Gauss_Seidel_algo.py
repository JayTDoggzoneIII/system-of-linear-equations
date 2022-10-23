from Matrix_Class import Matrix
from math import fabs
from sqrt_algo import sqrt
from find_eigenvalues_algo import find_eigenvalues

EPS = 1e-3
'''
def converge(A,x,b) -> bool:
    n = len(b)
    norm = 0
    for i in range(n):
        norm += (sum(A[i][k]*x[k] for k in range(n)) - b[i])**2
    return sqrt(norm) < EPS'''

converge_d = lambda A,x,b: (sum(fabs(sum(A[i][j]*x[j] for j in range(len(b))) - b[i]) for i in range(len(b)))) < EPS;

def converge(A:Matrix, x:list, b:list) -> bool:
    y = [sum(A[i][j]*x[j] for j in range(len(b))) for i in range(len(b))]
    for i in range(len(b)):
        if (fabs(y[i] - b[i]) > EPS): return False;
    return True;

def diag(A:Matrix) -> bool:
    n = len(A)
    for i in range(n):
        s = sum(fabs(A[i][j]) for j in range(n))
        s -= fabs(A[i][i])
        if (fabs(s - fabs(A[i][i])) > EPS): return False;
    return True;


def Gauss_Seidel(A:Matrix, b:list) -> list:
    n = len(A)
    A0 = Matrix([[A[i][j] for j in range(n)] for i in range(n)])
    b0 = [b[i] for i in range(n)]
    eig = find_eigenvalues(A)
    flag = True
    for i in range(n):
        if (type(eig[i]) == complex or eig[i] < -EPS):
            flag = False
            break;
    if (not flag):
        T = A.T()
        A = T*A
        b = [sum(T[i][k]*b[k] for k in range(n)) for i in range(n)] 
    
    C = Matrix(n)
    d = [0]*n
    for i in range(n):
        d[i] = b[i]/A[i][i]
        for j in range(n):
            C[i][j] = -A[i][j]/A[i][i] if (i != j) else 0
    x_prev = d[:]
    
    x = [sum(C[i][j]*x_prev[j] for j in range(i,n)) + d[i] for i in range(n)]
    
    count = 0
    while (not converge(A0,x,b0) and count <= 10**5):
        x_prev = x[:]
        
        for i in range(n):
            x[i] = sum(C[i][j]*x[j] for j in range(i)) + sum(C[i][j]*x_prev[j] for j in range(i,n)) + d[i]        
        
        count += 1
    return x,count;
