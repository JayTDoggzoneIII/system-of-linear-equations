from Matrix_Class import Matrix
from math import fabs
from sqrt_algo import sqrt
from Householder_algo import Householder_decomposition

def quasi_triangular(A:Matrix, EPS1 = 1e-9) -> bool:
    n = len(A)
    for i in range(n):
        for j in range(i-1):
            if (fabs(A[i][j]) > EPS1): return False;
            
    for i in range(n-1):
        if (fabs(A[i+1][i]) < EPS1): continue;
        else:
            flag1 = flag2 = True
            if (i + 2 < n): flag1 = fabs(A[i+2][i+1]) < EPS1
            if (0 < i - 1): flag2 = fabs(A[i][i-1]) < EPS1
            if (not (flag1 and flag2)): return False;
    return True;
'''
def some_check(A:Matrix, EPS1 = 1e-9) -> bool:
    n = len(A)
    k = 0
    for i in range(n):
        for j in range(i-1):
            k += fabs(A[i][j])
    return k < EPS1'''

def find_eigenvalues(A:Matrix) -> list:
    n = len(A)
    curA = Matrix([[A[i][j] for j in range(n)] for i in range(n)])
    nextA = Matrix([[A[i][j] for j in range(n)] for i in range(n)])
    count = 0
    while (True):
        if (quasi_triangular(nextA) or count > 50_000): break;
        curA = Matrix([[nextA[i][j] for j in range(n)] for i in range(n)])
        shift = Matrix([[curA[n-1][n-1]*int(i==j) for j in range(n)] for i in range(n)])
        Q,R = Householder_decomposition(curA - shift)
        nextA = R*Q + shift
        count += 1
    ans = [0]*n
    i = 0
    while (i < n):
        if (i + 1 == n or fabs(nextA[i+1][i]) < 1e-9): 
            ans[i] = nextA[i][i]
        else:
            Re_z = (nextA[i][i] + nextA[i+1][i+1])/2
            det = (nextA[i][i]*nextA[i+1][i+1] - nextA[i][i+1]*nextA[i+1][i])
            Im_z = (det - Re_z**2)**0.5
            ans[i],ans[i+1] = (Re_z + 1j*Im_z), (Re_z - 1j*Im_z)
            if (fabs(ans[i].imag) < 1e-11 and fabs(ans[i+1].imag) < 1e-11): ans[i],ans[i+1] = ans[i].real,ans[i+1].real
            i += 1
        i += 1
    return ans;
    #return nextA,curA,ans,count

