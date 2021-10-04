module phot
use iso_c_binding
use ellip

implicit none

real*8, parameter :: pi = 3.14159265358979323846, pihalf = 1.5707963267948966, twopithree = 2.0943951023931953
real*8, parameter :: o3 = 0.33333333333333333333, o9 = 0.1111111111111111111, twopi = 6.283185307179586
real*8, parameter :: pithird = 1.0471975511965976, pisixth = 0.5235987755982988

contains

subroutine compute_phis(rp, rm, bp2, bm2, bpm2, bp, bm, bpm, cosphi2, cosphi1, phi1, phi2)

    real*8 :: rp, rm, bp2, bm2, bpm2
    real*8 :: rp2, rm2, denom
    real*8, intent(out) :: cosphi1, cosphi2, phi1, phi2
    real*8 :: delta, gamma, s, t, tmp
    real*8 :: bp, bm, bpm
    real*8 :: a, b, c, theta1, theta2
    
    a = rm
    b = bpm
    c = rp
    if (a .gt. b) then
        a = bpm
        b = rm
    end if
    if (b .gt. c) then
        theta1 = 2 * Atan2(Sqrt(((a - b) + c) * (c - (a - b))), Sqrt((a + (b + c)) * ((a - c) + b)))
    else
        theta1 = 2 * Atan2(Sqrt(((a - b) + c) * (b - (a - c))), Sqrt((a + (b + c)) * ((a - c) + b)))
    end if
    
    
    a = bm
    b = bpm
    c = bp
    if (a .gt. b) then
        a = bpm
        b = bm
    end if
    if (b .gt. c) then
        theta2 = 2 * Atan(Sqrt(((a - b) + c) * (c - (a - b))), Sqrt((a + (b + c)) * ((a - c) + b)))
    else
        theta2 = 2 * Atan(Sqrt(((a - b) + c) * (b - (a - c))), Sqrt((a + (b + c)) * ((a - c) + b)))
    end if

    phi1 = theta2 + theta1
    phi2 = theta2 - theta1
    if (phi1 .gt. pi) then
        phi1 = phi1 - twopi
    end if
    if (phi2 .gt. pi) then
        phi2 = phi2 - twopi
    end if
    cosphi2 = Cos(phi2)
    cosphi1 = Cos(phi1)
end

! angle is with respect to the r1 circle's center
subroutine compute_theta(r1, r12, r2, r22, b1, b2, costheta, theta)

    real*8 :: a, b, c, r1, r2, r12, r22, b1, b2, costheta, theta

    a = b1
    b = r1
    c = r2
    if (a .gt. b) then
        a = r1
        b = b1
    end if
    if (b .gt. c) then
        theta = 2 * Atan(Sqrt(((a - b) + c) * (c - (a - b))), Sqrt((a + (b + c)) * ((a - c) + b)))
    else
        theta = 2 * Atan(Sqrt(((a - b) + c) * (b - (a - c))), Sqrt((a + (b + c)) * ((a - c) + b)))
    end if

    costheta = Cos(theta)
end

subroutine flux(c1, c2, rp, rm, bp2, bm2, bpm2, lc, j) bind(C, name="flux")

    integer (c_int), bind(C) :: j
    integer :: i
    real (c_double), bind(C) :: rp, rm
    real (c_double), bind(C), dimension(j) :: bp2, bm2, bpm2
    real (c_double), bind(C), intent(out), dimension(j) :: lc
    real (c_double), bind(C) :: c1, c2
    
    real*8 :: costhetapm, costhetapstar, costhetamstar, cosphim1, cosphim2, cosphip1, cosphip2, costheta, theta
    real*8 :: thetapstar, thetamstar, thetapm, phim1, phim2, phip1, phip2, thetap, thetam, costhetap, costhetam
    real*8 :: phim, phip, cosphim, cosphip, cosphi1, cosphi2, tmp
    real*8 :: bp2i, bm2i, bpm2i, bpi, bmi, bpmi, rp2, rm2, obpi, obmi
    real*8 :: d1, d2, a
    real*8 :: f0, of0
    
    real*8 :: bp(j), bm(j), bpm(j)
    bp = Sqrt(bp2)
    bm = Sqrt(bm2)
    bpm = Sqrt(bpm2)
    
    rp2 = rp * rp
    rm2 = rm * rm
    
    f0 = (1.d0 - c1 - 2 * c2) * pi + (c1 + 2 * c2) * twopithree + c2 * pihalf
    of0 = 1.d0 / f0
    lc = 1.d0
    
    do i=1,j,1
    
        bp2i = bp2(i)
        bm2i = bm2(i)
        bpm2i = bpm2(i)
        bpi = bp(i)
        bmi = bm(i)
        bpmi = bpm(i)
        
        if ((bpi .gt. rp + 1.d0) .AND. (bmi .gt. rm + 1.d0)) then
            lc(i) = 1.d0
            goto 1
        
        else if (bpmi .gt. rp + rm) then
            if (bpi .gt. rp + 1.d0) then
                if (bmi .gt. rm + 1.d0) then
                    lc(i) = 1.d0
                    goto 1
                else
                    if (bmi + rm .lt. 1.d0) then
                        lc(i) = 1.d0 - Arc(c1, c2, -pi, pi, 1.d0, 1.d0, rm, bmi) * of0
                        goto 1
                    else
                        call compute_theta(rm, rm2, 1.d0, 1.d0, bmi, bm2i, costhetamstar, thetamstar)
                        call compute_theta(1.d0, 1.d0, rm, rm2, bmi, bm2i, costheta, theta)
                        lc(i) = 2 * (F(c1, c2, pi - theta, -costheta, 1.d0, 0.d0) &
                                - F(c1, c2, thetamstar, costhetamstar, rm, bmi)) * of0
                        goto 1
                    end if
                end if
            else
                if (bmi .gt. rm + 1.d0) then
                    if (bpi + rp .lt. 1.d0) then
                        lc(i) = 1.d0 - 2 * F(c1, c2, pi, 1.d0, rp, bpi) * of0
                        goto 1
                    else
                        call compute_theta(rp, rp2, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                        call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costheta, theta)
                        lc(i) = 2 * (F(c1, c2, pi - theta, -costheta, 1.d0, 0.d0) &
                              - F(c1, c2, thetapstar, costhetapstar, rp, bpi)) * of0
                        goto 1
                    end if
                else
                    if (bpi + rp .lt. 1.d0) then
                        if (bmi + rm .lt. 1.d0) then
                            lc(i) = 1.d0 - 2 * (F(c1, c2, pi, 1.d0, rm, bmi) &
                                  + F(c1, c2, pi, 1.d0, rp, bpi)) * of0
                            goto 1
                        else
                            call compute_theta(rm, rm2, 1.d0, 1.d0, bmi, bm2i, costhetamstar, thetamstar)
                            call compute_theta(1.d0, 1.d0, rm, rm2, bmi, bm2i, costheta, theta)
                            lc(i) = 2 * (F(c1, c2, pi - theta, -costheta, 1.d0, 0.d0) &
                                  - F(c1, c2, thetamstar, costhetamstar, rm, bmi) &
                                  - F(c1, c2, pi, 1.d0, rp, bpi)) * of0
                            goto 1
                        end if
                    else
                        if (bmi + rm .lt. 1.d0) then
                            call compute_theta(rp, rp2, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                            call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costheta, theta)
                            lc(i) = 2 * (F(c1, c2, pi - theta, -costheta, 1.d0, 0.d0) &
                                  - F(c1, c2, thetapstar, costhetapstar, rp, bpi) & 
                                  - F(c1, c2, pi, 1.d0, rm, bmi)) * of0
                            goto 1
                        else
                            call compute_theta(rp, rp2, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                            call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costhetap, thetap)
                            call compute_theta(rm, rm2, 1.d0, 1.d0, bmi, bm2i, costhetamstar, thetamstar)
                            call compute_theta(1.d0, 1.d0, rm, rm2, bmi, bm2i, costhetam, thetam)
                            theta = thetap + thetam
                            lc(i) = 2 * (F(c1, c2, pi - theta, -Cos(theta), 1.d0, 0.d0) &
                                  - F(c1, c2, thetapstar, costhetapstar, rp, bpi) & 
                                  - F(c1, c2, thetamstar, costhetamstar, rm, bmi)) * of0
                            goto 1
                        end if
                    end if
                end if
            end if
        else
            if (bpi .gt. rp + 1.d0) then
                if (bmi .gt. rm + 1.d0) then
                    lc(i) = 1.d0
                    goto 1
                else
                    if (bmi + rm .lt. 1.d0) then
                        lc(i) = 1.d0 - 2 * F(c1, c2, pi, 1.d0, rm, bmi) * of0
                        goto 1
                    else
                        call compute_theta(rm, rm2, 1.d0, 1.d0, bmi, bm2i, costhetamstar, thetamstar)
                        call compute_theta(1.d0, 1.d0, rm, rm2, bmi, bm2i, costheta, theta)
                        lc(i) = 2 * (F(c1, c2, pi - theta, -costheta, 1.d0, 0.d0) &
                              - F(c1, c2, thetamstar, costhetamstar, rm, bmi)) * of0
                        goto 1
                    end if
                end if
            else
                if (bmi .gt. rm + 1.d0) then
                    if (bpi + rp .lt. 1.d0) then
                        lc(i) = 1.d0 - 2 * F(c1, c2, pi, 1.d0, rp, bpi) * of0
                        goto 1
                    else
                        call compute_theta(rp, rp2, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                        call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costheta, theta)
                        lc(i) = 2 * (F(c1, c2, pi - theta, -costheta, 1.d0, 0.d0) &
                              - F(c1, c2, thetapstar, costhetapstar, rp, bpi)) * of0
                        goto 1
                    end if
                else
                    if (bpi + rp .lt. 1.d0) then
                        if (bmi + rm .lt. 1.d0) then
                            if (bpmi + rm .lt. rp) then
                                lc(i) = 1.d0 - 2 * F(c1, c2, pi, 1.d0, rp, bpi) * of0
                                goto 1
                            else
                                call compute_phis(rp, rm, bp2i, bm2i, bpm2i, &
                                    bpi, bmi, bpmi, cosphim2, cosphim1, phim1, phim2)
                                call compute_phis(rm, rp, bm2i, bp2i, bpm2i, &
                                    bmi, bpi, bpmi, cosphip2, cosphip1, phip1, phip2)
                                lc(i) = 1.d0 - (Arc(c1, c2, phip1, phip2, cosphip1, cosphip2, rp, bpi) &
                                    + Arc(c1, c2, phim1, phim2, cosphim1, cosphim2, rm, bmi)) * of0
                                goto 1
                            end if
                        else 
                            if (bpmi + rm .lt. rp) then
                                goto 1
                            else
                                call compute_phis(rp, rm, bp2i, bm2i, bpm2i, &
                                    bpi, bmi, bpmi, cosphim2, cosphim1, phim1, phim2)
                                call compute_phis(rm, rp, bm2i, bp2i, bpm2i, &
                                    bmi, bpi, bpmi, cosphip2, cosphip1, phip1, phip2)
                                call compute_theta(rm, rm2, 1.d0, 1.d0, bmi, bm2i, costhetamstar, thetamstar)
                                call compute_theta(1.d0, 1.d0, rm, rm2, bmi, bm2i, costhetam, thetam)
                                lc(i) = (2 * F(c1, c2, pi - thetam, -costhetam, 1.d0, 0.d0) &
                                  - Arc(c1, c2, -thetamstar, phim2, costhetamstar, cosphim2, rm, bmi) &
                                  - Arc(c1, c2, phim1, thetamstar, cosphim1, costhetamstar, rm, bmi) &
                                  - Arc(c1, c2, phip1, phip2, cosphip1, cosphip2, rp, bpi)) * of0
                                goto 1
                            end if
                        end if
                    else
                        if (bmi + rm .lt. 1.d0) then 
                            if (bpmi + rm .lt. rp) then
                                call compute_theta(rp, rp2, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                                call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costheta, theta)
                                lc(i) = 2 * (F(c1, c2, pi - theta, -costheta, 1.d0, 0.d0) &
                                    - F(c1, c2, thetapstar, costhetapstar, rp, bpi)) * of0
                                goto 1
                            else
                                call compute_phis(rp, rm, bp2i, bm2i, bpm2i,&
                                    bpi, bmi, bpmi, cosphim2, cosphim1, phim1, phim2)
                                call compute_phis(rm, rp, bm2i, bp2i, bpm2i, &
                                    bmi, bpi, bpmi, cosphip2, cosphip1, phip1, phip2)
                                call compute_theta(rp, rp2, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                                call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costhetap, thetap)
                                lc(i) = (2 * F(c1, c2, pi - thetap, -costhetap, 1.d0, 0.d0) &
                                  - Arc(c1, c2, -thetapstar, phip2, costhetapstar, cosphip2, rp, bpi) &
                                  - Arc(c1, c2, phip1, thetapstar, cosphip1, costhetapstar, rp, bpi) &
                                  - Arc(c1, c2, phim1, phim2, cosphim1, cosphim2, rm, bmi)) * of0
                                goto 1
                            end if
                        else
                            obpi = 1.d0 / bpi
                            obmi = 1.d0 / bmi
                            if (bpmi + rm .lt. rp) then
                                call compute_theta(rp, rp2, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                                call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costheta, theta)
                                lc(i) = 2 * (F(c1, c2, pi - theta, -costheta, 1.d0, 0.d0) &
                                    - F(c1, c2, thetapstar, costhetapstar, rp, bpi)) * of0
                                goto 1
                            else
                                costhetapm = (bp2i + bm2i - bpm2i) * 0.5 * obpi * obmi
                                cosphim = (bm2i + 1.d0 - rm2) * 0.5 * obmi
                                cosphip = (bp2i + 1.d0 - rp2) * 0.5 * obpi
                                thetapm = Acos(costhetapm)
                                phim = Acos(cosphim)
                                phip = Acos(cosphip)
                                if (thetapm + phim .lt. phip) then
                                        call compute_phis(rp, rm, bp2i, bm2i, bpm2i, &
                                            bpi, bmi, bpmi, cosphim2, cosphim1, phim1, phim2)
                                        call compute_phis(rm, rp, bm2i, bp2i, bpm2i, &
                                            bmi, bpi, bpmi, cosphip2, cosphip1, phip1, phip2)
                                        call compute_theta(rp, rp2, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                                        call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costhetap, thetap)
                                        if (phip2 .gt. thetapstar) then
                                            lc(i) = 2 * (F(c1, c2, pi - thetap, -costhetap, 1.d0, 0.d0) &
                                                - F(c1, c2, thetapstar, costhetapstar, rp, bpi)) * of0
                                            goto 1
                                        else
                                            lc(i) = (2 * F(c1, c2, pi - thetap, -costhetap, 1.d0, 0.d0) &
                                              - Arc(c1, c2, -thetapstar, phip2, costhetapstar, cosphip2, rp, bpi) &
                                              - Arc(c1, c2, phip1, thetapstar, cosphip1, costhetapstar, rp, bpi) &
                                              - Arc(c1, c2, phim1, phim2, cosphim1, cosphim2, rm, bmi)) * of0
                                             goto 1
                                        end if
                                else if (thetapm + phip .lt. phim) then
                                    if ((bpi - rp) .lt. (bmi - rm)) then
                                        call compute_phis(rp, rm, bp2i, bm2i, bpm2i, &
                                            bpi, bmi, bpmi, cosphim2, cosphim1, phim1, phim2)
                                        call compute_phis(rm, rp, bm2i, bp2i, bpm2i, &
                                            bmi, bpi, bpmi, cosphip2, cosphip1, phip1, phip2)
                                        call compute_theta(rm, rm2, 1.d0, 1.d0, bmi, bm2i, costhetamstar, thetamstar)
                                        call compute_theta(1.d0, 1.d0, rm, rm2, bmi, bm2i, costhetam, thetam)
                                        lc(i) = (2 * F(c1, c2, pi - thetam, -costhetam, 1.d0, 0.d0) &
                                            - Arc(c1, c2, -thetamstar, phim2, costhetamstar, cosphim2, rm, bmi) &
                                            - Arc(c1, c2, phim1, thetamstar, cosphim1, costhetamstar, rm, bmi) & 
                                            - Arc(c1, c2, phip1, phip2, cosphip1, cosphip2, rp, bpi)) * of0
                                        goto 1
                                    else
                                        call compute_theta(rm, rm2, 1.d0, 1.d0, bmi, bm2i, costhetamstar, thetamstar)
                                        call compute_theta(1.d0, 1.d0, rm, rm2, bmi, bm2i, costheta, theta)
                                        lc(i) = 2 * (F(c1, c2, pi - theta, -costheta, 1.d0, 0.d0) &
                                          - F(c1, c2, thetamstar, costhetamstar, rm, bmi)) * of0
                                        goto 1
                                    end if
                                else
                                    costheta = (bpm2i + bm2i - bp2i) / (2 * bpmi * bmi)
                                    cosphim = (bpm2i + rm2 - rp2) / (2 * bpmi * rm)
                                    cosphi1 = Cos(Acos(costheta) - Acos(cosphim))
                                    cosphi2 = Cos(Acos(costheta) + Acos(cosphim))
                                    d1 = rm2 + bm2i - 2 * rm * bmi * cosphi1
                                    d2 = rm2 + bm2i - 2 * rm * bmi * cosphi2
                                    if (d1 .gt. 1.d0) then
                                        call compute_theta(rp, rp2, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                                        call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costhetap, thetap)
                                        call compute_theta(rm, rm2, 1.d0, 1.d0, bmi, bm2i, costhetamstar, thetamstar)
                                        call compute_theta(1.d0, 1.d0, rm, rm2, bmi, bm2i, costhetam, thetam)
                                        theta = thetap + thetam
                                        lc(i) = 2 * (F(c1, c2, pi - theta, -Cos(theta), 1.d0, 0.d0) &
                                          - F(c1, c2, thetapstar, costhetapstar, rp, bpi) & 
                                          - F(c1, c2, thetamstar, costhetamstar, rm, bmi)) * of0
                                        goto 1
                                    else if (d2 .lt. 1.d0) then
                                        call compute_phis(rp, rm, bp2i, bm2i, bpm2i, &
                                            bpi, bmi, bpmi, cosphim2, cosphim1, phim1, phim2)
                                        call compute_phis(rm, rp, bm2i, bp2i, bpm2i, &
                                            bmi, bpi, bpmi, cosphip2, cosphip1, phip1, phip2)
                                        call compute_theta(rm, rm2, 1.d0, 1.d0, bmi, bm2i, costhetamstar, thetamstar)
                                        call compute_theta(1.d0, 1.d0, rm, rm2, bmi, bm2i, costhetam, thetam)
                                        call compute_theta(rp, rp2, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                                        call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costhetap, thetap)
                                        theta = thetam + thetap
                                        lc(i) = (2 * F(c1, c2, pi - theta, -Cos(theta), 1.d0, 0.d0) &
                                            - Arc(c1, c2, -thetamstar, -phim1, costhetamstar, cosphim1, rm, bmi) & 
                                            - Arc(c1, c2, -phim2, thetamstar, cosphim2, costhetamstar, rm, bmi) &
                                            - Arc(c1, c2, phip1, thetapstar, cosphip1, costhetapstar, rp, bpi) &
                                            - Arc(c1, c2, -thetapstar, phip2, costhetapstar, cosphip2, rp, bpi)) * of0
                                        goto 1
                                    else
                                        call compute_phis(rp, rm, bp2i, bm2i, bpm2i, &
                                            bpi, bmi, bpmi, cosphim2, cosphim1, phim1, phim2)
                                        call compute_phis(rm, rp, bm2i, bp2i, bpm2i, &
                                            bmi, bpi, bpmi, cosphip2, cosphip1, phip1, phip2)
                                        call compute_theta(rm, rm2, 1.d0, 1.d0, bmi, bm2i, costhetamstar, thetamstar)
                                        call compute_theta(1.d0, 1.d0, rm, rm2, bmi, bm2i, costhetam, thetam)
                                        call compute_theta(rp, rp2, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                                        call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costhetap, thetap)
                                        theta = 0.5 * (thetap + thetam + thetapm)
                                        lc(i) = (2 * F(c1,c2, pi - theta, -Cos(theta), 1.d0, 0.d0) & 
                                            - Arc(c1, c2, -phim2, thetamstar, cosphim2, costhetamstar, rm, bmi) & 
                                            - Arc(c1, c2, -thetapstar, phip2, costhetapstar, cosphip2, rp, bpi)) * of0
                                        goto 1
                                    end if
                                end if
                            end if
                        end if
                    end if
                end if
            end if  
1        end if
    end do
    return
    
end

real*8 function Arc(c1, c2, phi1, phi2, cphi1, cphi2, r, b)

    real*8 :: phi1, phi2, r, b, c1, c2, cphi1, cphi2
    real*8 :: const, lin, quad
        
    if (phi1 < 0) then
        if (phi2 < 0) then
            if (phi2 < phi1) then
                Arc = 2 * F(c1, c2, pi, 1.d0, r, b) + F(c1, c2, -phi1, cphi1, r, b) - F(c1, c2, -phi2, cphi2, r, b)
                return
            else
                Arc = -F(c1, c2, -phi2, cphi2, r, b) + F(c1, c2, -phi1, cphi1, r, b)
                return
            end if
        else
            Arc = F(c1, c2, phi2, cphi2, r, b) + F(c1, c2, -phi1, cphi1, r, b)
            return
        end if
    else
        if (phi2 < 0) then
            Arc = 2 * F(c1, c2, pi, 1.d0, r, b) - F(c1, c2, phi1, cphi1, r, b) - F(c1, c2, -phi2, cphi2, r, b)
            return
        else
            if (phi2 < phi1) then
                Arc = 2 * F(c1, c2, pi, 1.d0, r, b) + F(c1, c2, phi2, cphi2, r, b) - F(c1, c2, phi1, cphi1, r, b)
                return
            else
                Arc = F(c1, c2, phi2, cphi2, r, b) - F(c1, c2, phi1, cphi1, r, b)
                return
            end if
        end if
    end if
    
    return
    
end function

real*8 function F(c1, c2, phi, cphi, r, b)

    real*8 :: c1, c2, phi, cphi, r, b, Fc, Fq, Fl
    real*8 :: o, ac, bc
    real*8 :: alpha, beta, gamma, d, s, n, m, x, y, ome, tans, sphi, br, bmr, bpr
    real*8 :: r2, b2, bmr2, tanphi2, sinphi2, apb, apbo
    real*8 :: ellippi, ellipf_tmp, ome_tmp, tans_tmp, eplusf
    
    r2 = r * r
    b2 = b * b
    bmr = b - r
    bmr2 = bmr * bmr
    bpr = b + r
    br = b * r
        
    if (phi .eq. pi) then
        Fc = r2 * pihalf     
        Fq = pihalf * r2 * (b2 + 0.5 * r2)        
        if (-c1 .eq. 2 * c2) then
            Fl = 0.d0
        else
            Fl = 0.d0
            if (b == 0.d0) then
                if (r == 1.d0) then
                    Fl = pithird
                    goto 2
                else
                    Fl = pithird * (1.d0 - (1.d0 - r2) ** (1.5))
                    goto 2
                end if
            else if (b == r) then
                if (r == 0.5) then
                    Fl = pisixth + 2.d0 * o9
                    goto 2
                else
                    alpha = 4 * (2 * r2 - 1.d0) * o9
                    beta = (1.d0 - 4 * r2) * o9
                    d = pisixth
                    s = pihalf
                    m = 4 * r2
                    ome = Sqrt(1.d0 - m)
                    eplusf = el2(Tan(s), ome, alpha + beta, beta + alpha * ome * ome)
                    Fl = eplusf + d
                    goto 2
                end if
            else if (bpr == 1.d0) then
                y = Sqrt(br)
                Fl = pisixth + Atan(2 * y / (1.d0 - 2 * r)) * o3 &
                    + pisixth *  Sqrt(1.d0 - 4 * r * b) / (2 * r - 1.d0) &
                    + 2.d0 * o9 * y * (2 * r * (5 * r - 2.d0) &
                    + 2 * br - 3.d0)
                goto 2
            else if (bpr .gt. 1.d0) then
        
                x = Sqrt(1.d0 - bmr2)
                y = Sqrt(br)
                alpha = 2 * y * (7 * r2 + b2 - 4.d0) * o9
                beta = 3.d0 + 2*r*(b2 * b + 5 * b2 * r + 3*r*(-2.d0 + r2) + b*(-4.d0 + 7*r2))
                beta = -beta * o9 * 0.5 / y
                gamma = (bpr) * o3 / (2 * (bmr) * y)
                m = (1.d0 - bmr2) / (4 * r * b)
                ome = Sqrt(1.d0 - m)
                if (abs(1.d0 / Sqrt(m) - 1.D0) .lt. 1.D-12) then
                    d = pisixth
                    n = 1.d0 - 1.d0 / bmr2
                    ac = 1.d0
                    bc = 1.d0
                    ome_tmp = ome
                    ellippi = cel(ome_tmp, 1.d0 - n, ac, bc)
                    o = 1.d0
                    ome_tmp = ome
                    eplusf = cel(ome_tmp, o, alpha + beta, beta + alpha * ome * ome)
                else
                    d = pisixth
                    tans = 1.d0 / Sqrt(m - 1.d0)
                    n = 1.d0 - 1.d0 / bmr2
                    tans_tmp = tans
                    ome_tmp = ome
                    ellippi = el3(tans_tmp, ome_tmp, 1.d0 - n)
                    tans_tmp = tans
                    eplusf = el2(tans_tmp, ome_tmp, alpha + beta, beta + alpha * ome * ome)
                end if
                Fl = eplusf + gamma * ellippi + d
                goto 2
            else
                x = Sqrt(1.d0 - bmr2)
                alpha = (7 * r2 + b2 - 4.d0) * x * o9
                beta = (r2*r2 + b2*b2 + r2 - b2 * (5.d0 + 2 * r2) + 1.d0) / (9.d0 * x)
                gamma = bpr / (bmr * 3.d0 * x)
                n = - 4 * r * b / bmr2
                m = 4 * r * b / (1.d0 - bmr2)
                ome = Sqrt(1.d0 - m)
                d = pisixth * (1.d0 - Sign(1.d0, bmr))
                ac = 1.d0
                bc = 1.d0
                ome_tmp = ome
                ellippi = cel(ome_tmp, 1.d0 - n, ac, bc) 
                o = 1.d0
                ome_tmp = ome
                eplusf = cel(ome_tmp, o, alpha + beta, beta + alpha * ome * ome)
                Fl = eplusf + gamma * ellippi + d
                goto 2
            end if
        end if
    else
        !sphi = Sqrt(1.d0 - cphi * cphi)
        !tanphi2 = sphi / (cphi + 1.d0)
        !sinphi2 = Sqrt((1.d0 - cphi) * 0.5)
        sphi = Sin(phi)
        sinphi2 = Sin(phi * 0.5)
        tanphi2 = Tan(phi * 0.5)
        
        Fc = 0.5 * (r2 * phi - br * sphi)
        Fq = 0.5 * (phi * r2 * (b2 + 0.5 * r2) - sphi * br * (b2 * o3 + 1.5 * r2)) &
            + 0.25 * sphi * r2 * (cphi * cphi * (br - r2 * cphi * o3) &
            + sphi * sphi * (r2 * cphi - br * o3))
        
        if (-c1 .eq. 2 * c2) then
            Fl = 0.d0
        else
            Fl = 0.d0
            if (b == 0.d0) then
                if (r == 1.d0) then
                    Fl = phi * o3
                    goto 2
                else
                    Fl = phi * (1.d0 - (1.d0 - r2) ** (1.5)) * o3
                    goto 2
                end if
            else if (b == r) then
                if (r == 0.5) then
                    Fl = phi * o3 * 0.5 + o3 * sinphi2 &
                       * (1.d0 - sinphi2 * sinphi2 * o3)
                    goto 2
                else
                    alpha = 4 * (2 * r2 - 1.d0) * o9
                    beta = (1.d0 - 4 * r2) * o9
                    d = phi * o3 * 0.5 - 2 * r2 * sphi & 
                        * Sqrt(1.d0 + 2 * r2 * (cphi-1.d0)) * o9
                    s = phi * 0.5
                    m = 4 * r2
                    ome = Sqrt(1.d0 - m)
                    eplusf = el2(Tan(s), ome, alpha + beta, beta + alpha * ome * ome)
                    Fl = eplusf + d
                    goto 2
                end if
            else if (bpr == 1.d0) then
                y = Sqrt(br)
                Fl = phi * o3 * 0.5 + Atan((2 * r - 1.d0) / tanphi2) * o3 &
                    + Atan(2 * y * sinphi2 / (1.d0 - 2 * r)) * o3 &
                    + pisixth *  Sqrt(1.d0 - 4 * r * b) / (2 * r - 1.d0) &
                    + 2.d0 * o9 * y * (2 * r * (5 * r - 2.d0) &
                    - 2 * br * cphi - 3.d0) * sinphi2
                goto 2
            else if (bpr .gt. 1.d0) then
        
                x = Sqrt(1.d0 - bmr2)
                y = Sqrt(br)
                alpha = 2 * y * (7 * r2 + b2 - 4.d0) * o9
                beta = 3.d0 + 2*r*(b2 * b + 5 * b2 * r + 3*r*(-2.d0 + r2) + b*(-4.d0 + 7*r2))
                beta = -beta * o9 * 0.5 / y
                gamma = (bpr) * o3 / (2 * (bmr) * y)
                m = (1.d0 - bmr2) / (4 * r * b)
                ome = Sqrt(1.d0 - m)
                if (abs(sinphi2 / Sqrt(m) - 1.D0) .lt. 1.D-12) then
                    d = phi * o3 * 0.5 - Atan((bpr / bmr) * tanphi2) * o3
                    n = 1.d0 - 1.d0 / bmr2
                    ac = 1.d0
                    bc = 1.d0
                    ome_tmp = ome
                    ellippi = cel(ome_tmp, 1.d0 - n, ac, bc)
                    o = 1.d0
                    ome_tmp = ome
                    eplusf = cel(ome_tmp, o, alpha + beta, beta + alpha * ome * ome)
                else
                    d = o3 * (phi * 0.5 - Atan((bpr / bmr) * tanphi2)) &
                        - (2 * br * o9) * sphi &
                        * Sqrt(1.d0 - b2 - r2 + 2 * br * cphi)
                    tans = 1.d0 / Sqrt(m / (sinphi2 * sinphi2) - 1.d0)
                    n = 1.d0 - 1.d0 / bmr2
                    tans_tmp = tans
                    ome_tmp = ome
                    ellippi = el3(tans_tmp, ome_tmp, 1.d0 - n)
                    tans_tmp = tans
                    eplusf = el2(tans_tmp, ome_tmp, alpha + beta, beta + alpha * ome * ome)
                end if
                Fl = eplusf + gamma * ellippi + d
                goto 2
            else
                x = Sqrt(1.d0 - bmr2)
                alpha = (7 * r2 + b2 - 4.d0) * x * o9
                beta = (r2*r2 + b2*b2 + r2 - b2 * (5.d0 + 2 * r2) + 1.d0) / (9.d0 * x)
                gamma = bpr / (bmr * 3.d0 * x)
                n = - 4 * r * b / bmr2
                m = 4 * r * b / (1.d0 - bmr2)
                ome = Sqrt(1.d0 - m)
                d = o3 * (phi * 0.5 - Atan((bpr / bmr) * tanphi2)) &
                    - 2 * br * o9 * sphi * Sqrt(1.d0 - b2 - r2 + 2 * br * cphi)
                ome_tmp = ome
                tans_tmp = tanphi2
                ellippi = el3(tans_tmp, ome_tmp, 1.d0 - n)
                tans_tmp = tanphi2
                eplusf = el2(tans_tmp, ome_tmp, alpha + beta, beta + alpha * ome * ome)
                Fl = eplusf + gamma * ellippi + d
                goto 2
            end if
        end if
    end if
    
2   F = (1.d0 - c1 - 2 * c2) * Fc + (c1 + 2 * c2) * Fl + c2 * Fq
    
    return
end function

end module phot