from math import fabs

def sqrt(n:float) -> float:
    x = 1
    while (True):
        nx = (x + n/x)/2
        if (fabs(x - nx) < 1e-14):
            break;
        x = nx;
    return x;