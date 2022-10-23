from Matrix_Class import Matrix
from Cube_solve import Kardano
from math import fabs
from sqrt_algo import sqrt
from find_eigenvalues_algo import find_eigenvalues

EPS = 1e-6

normF = lambda A: sqrt(sum(sum(A[i][j]*A[i][j] for j in range(len(A))) for i in range(len(A))));
norm1 = lambda A: max(sum(fabs(A[i][j]) for j in range(len(A))) for i in range(len(A)));


def spectral_norm(A:Matrix) -> float:
    eig = find_eigenvalues(A.T()*A)
    m = sqrt(eig[0].real**2 + eig[0].imag**2)
    for i in range(1,len(eig)):
        if (sqrt(eig[i].real**2 + eig[i].imag**2) > m): m = sqrt(eig[i].real**2 + eig[i].imag**2)
    return sqrt(m);

def diag(A:Matrix) -> bool:
    n = len(A)
    for i in range(n):
        s = sum(fabs(A[i][j]) for j in range(n))
        s -= fabs(A[i][i])
        if (s > fabs(A[i][i])): return False;
    return True;

def converge(sb,x,x_prev,sc,k) -> bool:
    norm = sqrt(sum((x[i] - x_prev[i])**2 for i in range(len(x))))
    return sb/(1-sb)*norm < EPS or (sb**k*sc)/(1 - sb) < EPS;

converge_t = lambda A,x,b: sqrt(sum((sum(A[i][j]*x[j] for j in range(len(b))) - b[i])**2 for i in range(len(b)))) < EPS;



def fixed_point_iteration(A:Matrix, b:list) -> list:
    n = len(A)
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
    k = 1/spectral_norm(A)
    c = [b[i]*k for i in range(n)]
    B = Matrix(n)
    for i in range(n):
        B[i][i] = 1
    B = B - A*k
    sb = spectral_norm(B)
    sc = sqrt(sum(c[i]*c[i] for i in range(len(c))))
    if (normF(B) >= 1 and norm1(B) >= 1 and sb >= 1):
        print("Botva")
        x_prev = c[:]
        x = [sum(B[i][j]*x_prev[j] for j in range(n)) + c[i] for i in range(n)]
        count = 0
        while (not converge_t(A,x_prev,b) and count <= 100_000):
            x_prev = x[:]
            x = [sum(B[i][j]*x_prev[j] for j in range(n)) + c[i] for i in range(n)]
            count += 1
        return x,count;
    x_prev = c[:]
    x = [sum(B[i][j]*x_prev[j] for j in range(n)) + c[i] for i in range(n)]
    count = 0
    while (not converge(sb,x,x_prev,sc,count) and count <= 100_000):

        x_prev = x[:]
        x = [sum(B[i][j]*x_prev[j] for j in range(n)) + c[i] for i in range(n)]
        count += 1
    return x,count;


