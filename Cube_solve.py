from decimal import Decimal as LD, getcontext
from cmath import acos, cos, pi
from math import inf
getcontext().prec = 50

def _cbrt(n:float) -> float:
    if (n<0):
        n = LD(str(n)[1:])
        return float(-(n**LD('0.'+'3'*50)));
    n = LD(str(n))
    return float(n**LD('0.'+'3'*50));    

def Kardano(a:float, b:float, c:float, d:float) -> tuple:
    if (not a):
        if (not b):
            if (not c):
                if (not d):
                    return +inf;
                return "No roots";
            return -d/c;
        D = c*c - 4*b*d
        if (D > 0):
            return (-c + D**0.5)/(2*b),(-c - D**0.5)/(2*b);
        if (not D):
            return -c/(2*b);
        return (-c + (1j * abs(D)**0.5))/(2*b),(-c - (1j * abs(D)**0.5))/(2*b);
    a,b,c = b/a,c/a,d/a
    p = -(a*a/3)+b
    q = 2*(a/3)*(a/3)*(a/3) - (a*b/3) + c
    Q = (p/3)*(p/3)*(p/3) + (q/2)*(q/2)
    if (not Q):
        A = _cbrt(-q/2 + Q**0.5)
        B = _cbrt(-q/2 - Q**0.5)        
        return 2*A-a/3,A-a/3,A-a/3;
    if (Q > 0):
        A = _cbrt(-q/2 + Q**0.5)
        B = _cbrt(-q/2 - Q**0.5)
        x1 = A+B - a/3
        x2 = -(A+B)/2 + (1j*((A-B)/2)*(3**0.5))-a/3
        if (x2.real < 1e-15):
            x2 = (0+x2.imag*1j)
        x3 = -(A+B)/2 - (1j*((A-B)/2)*(3**0.5))-a/3
        if (x3.real < 1e-15):
            x3 = (0+x3.imag*1j)       
        return round(x1,13),x2,x3;
    cosphi = -(q/2)*(3/-p)**(1.5)
    phi = acos(cosphi)
    return (2*(-p/3)**0.5 * cos(phi/3) - a/3).real,(2*(-p/3)**0.5 * cos(phi/3 + 2*pi/3) - a/3).real,(2*(-p/3)**0.5 * cos(phi/3 - 2*pi/3) - a/3).real;