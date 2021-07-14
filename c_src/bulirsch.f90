module bulirsch
use iso_c_binding
implicit none

contains

subroutine f_bulirsch(phi, m, e, j) bind(C, name="f_burl")

    integer (c_int), bind(C) :: j
    real (c_double), bind(C) :: phi(j), m(j)
    real (c_double), bind(C), dimension(j), intent(out) :: e
    
    call el1(tan(phi), sqrt(1 - m), e, j)
    return
end

subroutine e_bulirsch(phi, m, e, j) bind(C, name="e_burl")

    integer (c_int), bind(C) :: j
    real (c_double), bind(C) :: phi(j), m(j)
    real*8:: o(j)
    real (c_double), bind(C), dimension(j), intent(out) :: e
    
    o = 1.d0
    call el2(tan(phi), sqrt(1 - m), o, 1 - m, e, j)
    return
end

subroutine el1(x, kc, e, j) bind(C, name="el1")

    integer (c_int), bind(C) :: j
    real (c_double), bind(C) :: x(j), kc(j)
    real (c_double), bind(C), dimension(j), intent(out) :: e
    
    real*8 :: g(j), m(j), y(j)
    real*8 :: yi, ei, mi, kci, gi, xi
    real*8 :: D = 8
    real*8 :: ca, cb
    real*8 :: pi = 3.1415926535897932384626433832795
    integer :: l, i
    
    ca = 10**(-D/2.d0)
    cb = 10**(-D+2.d0)
    
    y = abs(1.d0/x)
    kc = abs(kc)
    m = 1.d0
    
    do i=1,j,1
    
        yi = y(i)
        ei = e(i)
        mi = m(i)
        kci = kc(i)
        gi = g(i)
        xi = x(i)
        l = 0
        
1       ei = mi * kci
        gi = mi
        mi = kci + mi
        yi = - (ei/yi) + yi
        if (yi == 0) then
            yi = sqrt(ei) * cb
        end if 
        if (abs(gi - kci) .gt. ca * gi) then
            kci = 2 * sqrt(ei)
            l = 2 * l
            if (yi .lt. 0) then
                l = l + 1
            end if 
            goto 1
        end if
        if (yi .lt. 0) then 
            l = l + 1
        end if 
        ei = (atan(mi/yi) + pi * l) / mi
        if (xi .lt. 0) then
            e(i) = -ei
        else
            e(i) = ei
        end if
    enddo
    
    return
end

subroutine el2(x, kc, a, b, e, j) bind(C, name="el2")

    integer (c_int), bind(C) :: j
    real (c_double), bind(C) :: x(j), kc(j), a(j), b(j)
    real (c_double), bind(C), dimension(j), intent(out) :: e
    
    real*8 :: c(j), dd(j), f(j), g(j), ik(j), m(j), p(j), y(j), z(j)
    real*8 :: ai, bi, ci, di, fi, iki, ppi, yi, zi, ei, mi, kci, gi, xi
    real*8 :: D = 8
    real*8 :: ca, cb
    real*8 :: pi = 3.1415926535897932384626433832795
    integer :: l, i
    
    ca = 10**(-D/2.d0)
    cb = 10**(-D+2.d0)
    
    c = x**2
    dd = c + 1.d0
    p = sqrt((1.d0 + kc**2 * c)/dd)
    dd = x / dd
    c = dd / (2.d0 * p)
    z = a - b
    ik = a
    a = (b + a) / 2.d0
    y = abs(1.d0/x)
    f = 0.d0
    kc = abs(kc)
    m = 1.d0
    
    do i=1,j,1
    
        ai = a(i)
        bi = b(i)
        ci = c(i)
        di = dd(i)
        ppi = p(i)
        zi = z(i)
        iki = ik(i)
        fi = f(i)
        yi = y(i)
        ei = e(i)
        mi = m(i)
        kci = kc(i)
        gi = g(i)
        xi = x(i)
        l = 0
        
1       bi = iki * kci + bi
        ei = mi * kci
        gi = ei / ppi
        di = fi * gi + di
        fi = ci
        iki = ai
        ppi = gi + ppi
        ci = (di / ppi + ci) / 2.d0
        gi = mi
        mi = kci + mi
        ai = (bi / mi + ai) / 2.d0
        yi = - (ei/yi) + yi
        if (yi == 0) then
            yi = sqrt(ei) * cb
        end if 
        if (abs(gi - kci) .gt. ca * gi) then
            kci = 2 * sqrt(ei)
            l = 2 * l
            if (yi .lt. 0) then
                l = l + 1
            end if 
            goto 1
        end if
        if (yi .lt. 0) then 
            l = l + 1
        end if 
        ei = (atan(mi/yi) + pi * l) * (ai / mi)
        if (xi .lt. 0) then
            e(i) = -ei
        else
            e(i) = ei + ci * zi
        end if
    enddo
    
    return
end

end module bulirsch