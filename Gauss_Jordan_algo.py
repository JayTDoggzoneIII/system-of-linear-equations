from Matrix_Class import Matrix
from math import fabs
INF = 10**9
EPS = 1e-11

def Gauss_Jordan(A:Matrix, b:list) -> list:
    n = len(b)
    
    a = [[0]*(n+1) for i in range(n)]
    for i in range(n):
        for j in range(n):
            a[i][j] = A[i][j]
        a[i][n] = b[i]
    
    where = [-1]*n
    col = row = 0
    while (col < n and row < n):
        sel = row
        for i in range(row,n):
            if (fabs(a[i][col]) > fabs(a[sel][col])): sel = i
        if (fabs(a[sel][col]) < EPS): continue;
        for i in range(col,n+1):
            a[sel][i],a[row][i] = a[row][i],a[sel][i]
        where[col] = row
        
        for i in range(n):
            if (i != row):
                c = a[i][col]/a[row][col]
                for j in range(col,n+1):
                    a[i][j] -= a[row][j]*c
        col += 1
        row += 1
    ans = [0]*n
    for i in range(n):
        if (where[i] != -1): ans[i] = a[where[i]][n] / a[where[i]][i]
    
 
    for i in range(n):
        s = sum(ans[j] * a[i][j] for j in range(n))
        if (fabs(s - a[i][n]) > EPS): return [None]*n;
        
    for i in range(n):
        if (where[i] == -1):
            return [INF]*n;
    return ans;