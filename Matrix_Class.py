from math import fabs


class Matrix:
    A = [[]]
    n = 0
    def __init__(self, x):
        if (type(x) == int):
            self.A = [[0]*x for i in range(x)]
            self.n = x
        elif (type(x) == list and len(x) and type(x[0]) == list):
            self.A = x[:]
            self.n = len(x)
    def __eq__(self,B) -> bool:
        for i in range(self.n):
            for j in range(self.n):
                if (fabs(self.A[i][j] - B[i][j]) > 1e-9): return False;
        return True;
            
    def __len__(self) -> int:
        return self.n;
    
    def __getitem__(self, i:int):
        return self.A[i];
    
    def __add__(self, B):
        if (self.n != len(B)): raise ValueError;
        ans = Matrix(self.n)
        for i in range(self.n):
            for j in range(self.n):
                ans[i][j] = self.A[i][j] + B[i][j]
        return ans;
    
    def __sub__(self, B):
        if (self.n != len(B)): raise ValueError;
        ans = Matrix(self.n)
        for i in range(self.n):
            for j in range(self.n):
                ans[i][j] = self.A[i][j] - B[i][j]
        return ans;
    
    def __mul__(self, B):
        ans = Matrix(self.n)
        if (type(B) == int or type(B) == float):
            for i in range(self.n):
                for j in range(self.n):
                    ans[i][j] = self.A[i][j]*B
            return ans;
        if (self.n != len(B)): raise ValueError;
        ans = Matrix(self.n)
        for i in range(self.n):
            for j in range(self.n):
                ans[i][j] = sum(self.A[i][k]*B[k][j] for k in range(self.n))
        return ans;
    
    def det(self) -> float:
        ans = 1
        C = Matrix(self.n)
        for i in range(self.n):
            for j in range(self.n):
                C[i][j] = self.A[i][j]
        for i in range(self.n):
            k = i
            for j in range(i+1, self.n):
                if (fabs(C[i][j]) > fabs(C[k][i])): k = j
                if (fabs(C[k][i]) < 1e-13): return 0;
            
            C.A[i], C.A[k] = C.A[k], C.A[i]
            
            if (i != k): ans = -ans
            ans *= C[i][i]
            
            for j in range(i+1, self.n): C[i][j] /= C[i][i]
            for j in range(self.n):
                if (i != j and fabs(C[i][i]) > 1e-9):
                    for k in range(i+1, self.n): C[j][k] -= C[i][k]*C[j][i]
        return ans;
    
    def T(self):
        return Matrix([[self.A[j][i] for j in range(self.n)] for i in range(self.n)]);

    def __repr__(self) -> str:
        return str(self.A).replace("], [", "]\n [");
    def __str__(self) -> str:
        return str(self.A);
    
    def inv(self):
        if (fabs(self.det()) < 1e-13): raise ValueError;
        C = [[0]*(2*self.n) for i in range(self.n)]
        for i in range(self.n):
            for j in range(self.n):
                C[i][j] = self.A[i][j]
                if (i == j): C[i][j+self.n] = 1
        for i in range(self.n):
            for j in range(self.n):
                if (i != j):
                    ratio = C[j][i]/C[i][i]
                    for k in range(2*self.n):
                        C[j][k] -= ratio*C[i][k]
        for i in range(self.n):
            div = C[i][i]
            for j in range(2*self.n):
                C[i][j] /= div
        return Matrix([[C[i][j] for j in range(self.n, 2*self.n)] for i in range(self.n)]);