module phot

use ellip
use iso_c_binding
implicit none

!PI=4.D0*DATAN(1.D0)
contains

!subroutine integrate_const(s1, s2, r, b, f) bind(C, name="integrate_const")
!    return
!end

subroutine integrate_linear(s1, s2, r, b, g) bind(C, name="integrate_linear")

    real (C_DOUBLE), bind(C) :: s1, s2, r, b
    real (C_DOUBLE), bind(C), intent(out) :: g
    real*8 :: d2, k2, k, c1, c2, c3, s3, t2, t1
    real*8 :: coeffA, ne, num, denom, coeffB
    real*8 :: R1, R2, DR, alpha, beta, gamma
    real*8 :: bint, dint, jint

    d2 = (r - b)**2
    k2 = 0.25 * (1.d0 - d2) / (r * b)
    ne = (d2 - 1) / d2
    
    c1 = sqrt(1.d0 - s1**2)
    c2 = sqrt(1.d0 - s2**2)
    
    c3 = c2 * c1 + s2 * s1 * sqrt((1.d0 - k2 * s2**2) & 
             * (1.d0 - k2 * s1**2)) / (1.d0 - k2 * s1**2 * s2**2)
    s3 = sqrt(1.d0 - c3**2)
    coeffA = - k2 * s1 * s2 * s3
    
    num = s1 * s2 * s3 * sqrt( ne * (1 - ne) * ( ne - k2))
    denom = 1 - ne * s3 + ne * s1 * s2 * s3 * sqrt(1 - k2 * s3)
    coeffB = sqrt( ne / ((1 - ne) * ( ne - k2 ))) * atan(num / denom)
    
    ! if something is amiss, check this logic. 
    if (abs(s1) .LT. 1e-10) then
        t1 = 0.d0
    else
        t1 = (1.d0 / 3) * atan(((r-b)/(r+b)) * (1 - c1) / s1)
    end if
    if (abs(s2) .LT. 1e-10) then
        t2 = 0.d0
    else
        t2 = (1.d0 / 3) * atan(((r-b)/(r+b)) * (1 - c2) / s2)
    end if
        
    if (b .GT. r) then
        R2 = -t2 + asin(s2)/6 + (2.d0 / 9) * r * b * sqrt(1 - r**2 - b**2 - 2 * r * b * c2) * s2
        R1 = -t1 + asin(s1)/6 + (2.d0 / 9) * r * b * sqrt(1 - r**2 - b**2 - 2 * r * b * c1) * s1
        DR = R2 - R1
    else if (b .LT. r) then
        R2 = t2 + asin(s2)/6 + (2.d0 / 9) * r * b * sqrt(1 - r**2 - b**2 - 2 * r * b * c2) * s2
        R1 = t1 + asin(s1)/6 + (2.d0 / 9) * r * b * sqrt(1 - r**2 - b**2 - 2 * r * b * c1) * s1
        DR = R2 - R1
    else 
        DR = 0
    end if
    
    alpha = (1.d0 / 6) * (1 - 4 * r**2 - 2 * r**4) / sqrt(r * b) &
        + sqrt(r * b) * (1.d0 / 9) * (7 * r**2 + 5 * r * b + b**2 - 4)
    beta = sqrt(r * b) * (1.d0 / 9) * (8 - 14 * r**2 - 2 * b**2)
    
    ! also check this limit... 
    if (r .EQ. b) then 
        gamma = 0
    else 
        gamma = (1.d0 / 6) * (r + b) / (r - b) / sqrt(r * b)
    end if 
    
    call elsbdj(s3, ne, 1.d0 - k2, bint, dint, jint)
    
    !f = (dR - beta * coeffA - gamma * coeffB) & 
    !    + (alpha + beta + gamma) * bint & 
    !    + (alpha + beta * (1 - k2) + gamma) * dint & 
    !    + gamma * ne * jint
    !f = jint   
    return
end

!subroutine integrate_quad(s1, s2, r, b, f) bind(C, name="integrate_quad")
!    return
!end 
    
end