module phot
use iso_c_binding
use ellip
implicit none

contains

real*8 function F(phi, r, b)

    real*8 :: phi, r, b
    real*8 :: pi = 3.14159265358979323846264
    real *8 :: o = 1.d0
    real*8 :: alpha, beta, gamma, d, s, n, m, x
    real*8 :: ellipf, ellipe, ellippi
    
    if (b == 0) then
        if (r == 0) then
            F = phi / 3.d0
            return
        else
            F = phi * (1.d0 - (1.d0 - r * r) ** (3.d0 * 0.5)) / 3.d0
            return
    else if (b == r) then
        if (r == 0.5) then
            F = phi / 6.d0 - (1 / 3.d0) * Sin(phi * 0.5) &
                   * (1.d0 - Sin(phi * 0.5)**2.d0 / 3.d0)
            return
        else
            alpha = 4 * (2 * r * r - 1.d0) / 9.d0
            beta = (1.d0 - 4 * r * r) / 9.d0
            d = phi / 6.d0 - 2 * r * r * Sin(phi) & 
                * Sqrt(1.d0 + 2 * r * r * (Cos(phi)-1.d0)) / 9.d0
            s = phi / 2.d0
            m = 4 * r * r
            ellipf = el1(Tan(s), Sqrt(1.d0 - m))
            ellipe = el2(Tan(s), Sqrt(1.d0 - m), o, 1.d0 - m)
            F = alpha * ellipf + beta * ellipe + d
            return
    else if (b + r == 1.d0) then
        F = phi / 6.d0 - Atan((2 * r - 1.d0) / Tan(phi * 0.5)) / 3.d0 &
                + Atan(2 * Sqrt(r * b) * Sin(phi * 0.5) / (1.d0 - 2 * r)) / 3.d0 &
                + pi *  Sqrt(1.d0 - 4 * r * b) / (2 * r - 1.d0) / 6.d0 &
                + (2.d0 / 9.d0) * Sqrt(b * r) * (2 * r * (5 * r - 2.d0) &
                - 2 * b * r * Cos(phi) - 3.d0) * Sin(phi * 2.d0)
        return
    else
        x = Sqrt(1.d0 - (b - r)**2.d0)
        alpha = (7 * r * r + b * b - 4.d0) * x / 9.d0
        beta = (r**4.d0 + b**4.d0 + r * r - b * b * (5.d0 + 2 * r * r) + 1.d0) / (9 * x)
        gamma = (b + r) / (b - r) / (3 * x)
        d = phi / 6.d0 - Atan((b + r) / (b - r) * Tan(phi * 0.5)) / 3.d0 &
                - (2 * b * r / 9.d0) * Sin(phi) &
                * Sqrt(1.d0 - b * b - r * r + 2 * b * r * Cos(phi))
        s = phi * 0.5
        n = - 4 * r * b / (b - r)**2
        m = 4 * r * b / (1.d0 - (r - b)**2)
        ellipf = el1(Tan(s), Sqrt(1.d0 - m))
        ellipe = el2(Tan(s), Sqrt(1.d0 - m), o, 1.d0 - m)
        ellippi = el3(Tan(s), sqrt(1.d0 - m), 1.d0 - n)
        F = alpha * elipf + beta * ellipe + gamma * ellippi + d
        return
        
    F = 0
    return

end module phot