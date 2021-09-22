module phot
use iso_c_binding
use ellip

implicit none

real*8, parameter :: pi = 3.14159265358979323846
real*8, parameter :: o3 = 0.33333333333333333333, o9 = 0.1111111111111111111

contains

subroutine compute_phis(rp, rm, bp2, bm2, bpm2, bp, bm, bpm, cosphi2, cosphi1, phi1, phi2)

    real*8 :: rp, rm, bp2, bm2, bpm2
    real*8 :: rp2, rm2
    real*8, intent(out) :: cosphi1, cosphi2, phi1, phi2
    real*8 :: delta, gamma, s, t, tmp
    real*8 :: bp, bm, bpm
    
    rp2 = rp * rp
    rm2 = rm * rm
    
    delta = bm * bpm2 * rm &
          * Sqrt(((bm - bp - bpm) * (bm + bp - bpm) * (bm - bp + bpm) &
          * (bm + bp + bpm) * (bpm - rm - rp) * (bpm + rm - rp) &
          * (bpm - rm + rp) * (bpm + rm + rp)) / (bm2 * bpm2 * bpm2 * rm2))
          
    gamma = bm2 * (bpm2 + rm2 - rp2) + (-bp2 + bpm2)*(bpm2 + rm2 - rp2)
    cosphi2 = (gamma - delta) / (4. * bm * bpm2 * rm)
    cosphi1 = (gamma + delta) / (4. * bm * bpm2 * rm)
    phi1 = Acos(cosphi1)
    phi2 = Acos(cosphi2)
    
    s = (bm2 * bm2 + (bp2 - bpm2)**2. & 
      - (bm2*(bpm2 * bpm2 + 2. * bp2 * rm2 & 
      + rm2 * rm2 - 2. * (bpm2 + rm2) * rp2 + rp2 * rp2)) / rm2) & 
      / (4. * bm2 * bpm2)
    s = Sign(1.d0, s)
    t = bpm * bpm - (-(bm2 * rm) + bp2 * rm - bm * rm2 + bm * rp2) & 
      / (bm + rm)
    t = Sign(1.d0, t)
    phi1 = s * t * phi1
    phi2 = - t * phi2
end

! angle is with respect to the r1 circle's center
subroutine compute_theta(r1, r12, r2, r22, b, b2, costheta, theta)

    real*8 :: r1, r2, b, costheta, theta, r12, r22, b2
    
    costheta = (r12 + b2 - r22) / (2.d0 * r1 * b)
    theta = Acos(costheta)
end

subroutine flux(c1, c2, rp, rm, bp2, bm2, bpm2, lc, j) bind(C, name="flux")

    integer (c_int), bind(C) :: j
    integer :: i
    real (c_double), bind(C) :: rp, rm
    real (c_double), bind(C) :: bp2(j), bm2(j), bpm2(j)
    real (c_double), bind(C), intent(out) :: lc(j)
    real (c_double), bind(C) :: c1, c2
    
    real*8 :: costhetapm, costhetapstar, costhetamstar, cosphim1, cosphim2, cosphip1, cosphip2, costheta, theta
    real*8 :: thetapstar, thetamstar, thetapm, phim1, phim2, phip1, phip2, thetap, thetam, costhetap, costhetam
    real*8 :: phim, phip, cosphim, cosphip, cosphi1, cosphi2
    real*8 :: bp2i, bm2i, bpm2i, bpi, bmi, bpmi, rp2, rm2
    real*8 :: d1, d2, a
    real*8 :: f0, of0
    
    real*8 :: bp(j), bm(j), bpm(j)
    bp = Sqrt(bp2)
    bm = Sqrt(bm2)
    bpm = Sqrt(bpm2)
    
    rp2 = rp * rp
    rm2 = rm * rm
    
    f0 = (1.d0 - c1 - 2 * c2) * pi + (c1 + 2 * c2) * 2 * pi * o3 + c2 * 0.5 * pi
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
                        lc(i) = 1.d0 - Arc(c1, c2, -pi, pi, rm, bmi) * of0
                        goto 1
                    else
                        call compute_theta(rm, rm2, 1.d0, 1.d0, bmi, bm2i, costhetamstar, thetamstar)
                        call compute_theta(1.d0, 1.d0, rm, rm2, bmi, bm2i, costheta, theta)
                        lc(i) = 2 * (F(c1, c2, pi - theta, 1.d0, 0.d0) &
                                - F(c1, c2, thetamstar, rm, bmi)) * of0
                        goto 1
                    end if
                end if
            else
                if (bmi .gt. rm + 1.d0) then
                    if (bpi + rp .lt. 1.d0) then
                        lc(i) = 1.d0 - 2 * F(c1, c2, pi, rp, bpi) * of0
                        goto 1
                    else
                        call compute_theta(rp, rp2, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                        call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costheta, theta)
                        lc(i) = 2 * (F(c1, c2, pi - theta, 1.d0, 0.d0) &
                              - F(c1, c2, thetapstar, rp, bpi)) * of0
                        goto 1
                    end if
                else
                    if (bpi + rp .lt. 1.d0) then
                        if (bmi + rm .lt. 1.d0) then
                            lc(i) = 1.d0 - 2 * (F(c1, c2, pi, rm, bmi) &
                                  + F(c1, c2, pi, rp, bpi)) * of0
                            goto 1
                        else
                            call compute_theta(rm, rm2, 1.d0, 1.d0, bmi, bm2i, costhetamstar, thetamstar)
                            call compute_theta(1.d0, 1.d0, rm, rm2, bmi, bm2i, costheta, theta)
                            lc(i) = 2 * (F(c1, c2, pi - theta, 1.d0, 0.d0) &
                                  - F(c1, c2, thetamstar, rm, bmi) &
                                  - F(c1, c2, pi, rp, bpi)) * of0
                            goto 1
                        end if
                    else
                        if (bmi + rm .lt. 1.d0) then
                            call compute_theta(rp, rp2, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                            call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costheta, theta)
                            lc(i) = 2 * (F(c1, c2, pi - theta, 1.d0, 0.d0) &
                                  - F(c1, c2, thetapstar, rp, bpi) & 
                                  - F(c1, c2, pi, rm, bmi)) * of0
                            goto 1
                        else
                            call compute_theta(rp, rp2, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                            call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costhetap, thetap)
                            call compute_theta(rm, rm2, 1.d0, 1.d0, bmi, bm2i, costhetamstar, thetamstar)
                            call compute_theta(1.d0, 1.d0, rm, rm2, bmi, bm2i, costhetam, thetam)
                            theta = thetap + thetam
                            lc(i) = 2 * (F(c1, c2, pi - theta, 1.d0, 0.d0) &
                                  - F(c1, c2, thetapstar, rp, bpi) & 
                                  - F(c1, c2, thetamstar, rm, bmi)) * of0
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
                        lc(i) = 1.d0 - 2 * F(c1, c2, pi, rm, bmi) * of0
                        goto 1
                    else
                        call compute_theta(rm, rm2, 1.d0, 1.d0, bmi, bm2i, costhetamstar, thetamstar)
                        call compute_theta(1.d0, 1.d0, rm, rm2, bmi, bm2i, costheta, theta)
                        lc(i) = 2 * (F(c1, c2, pi - theta, 1.d0, 0.d0) &
                              - F(c1, c2, thetamstar, rm, bmi)) * of0
                        goto 1
                    end if
                end if
            else
                if (bmi .gt. rm + 1.d0) then
                    if (bpi + rp .lt. 1.d0) then
                        lc(i) = 1.d0 - 2 * F(c1, c2, pi, rp, bpi) * of0
                        goto 1
                    else
                        call compute_theta(rp, rp2, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                        call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costheta, theta)
                        lc(i) = 2 * (F(c1, c2, pi - theta, 1.d0, 0.d0) &
                              - F(c1, c2, thetapstar, rp, bpi)) * of0
                        goto 1
                    end if
                else
                    if (bpi + rp .lt. 1.d0) then
                        if (bmi + rm .lt. 1.d0) then
                            if (bpmi + rm .lt. rp) then
                                lc(i) = 1.d0 - 2 * F(c1, c2, pi, rp, bpi) * of0
                                goto 1
                            else
                                call compute_phis(rp, rm, bp2i, bm2i, bpm2i, &
                                    bpi, bmi, bpmi, cosphim2, cosphim1, phim1, phim2)
                                call compute_phis(rm, rp, bm2i, bp2i, bpm2i, &
                                    bmi, bpi, bpmi, cosphip2, cosphip1, phip1, phip2)
                                lc(i) = 1.d0 - (Arc(c1, c2, phip1, phip2, rp, bpi) &
                                      + Arc(c1, c2, phim1, phim2, rm, bmi)) * of0
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
                                lc(i) = (2 * F(c1, c2, pi - thetam, 1.d0, 0.d0) &
                                      - Arc(c1, c2, -thetamstar, phim2, rm, bmi) &
                                      - Arc(c1, c2, phim1, thetamstar, rm, bmi) &
                                      - Arc(c1, c2, phip1, phip2, rp, bpi)) * of0
                                goto 1
                            end if
                        end if
                    end if
                                
            
                    if (bpi + rp .lt. 1.d0) then
                        if (bmi + rm .lt. 1.d0) then
                            if (bpmi + rm .lt. rp) then
                                lc(i) = 1.d0 - 2 * F(c1, c2, pi, rp, bpi) * of0
                                goto 1
                            else
                                call compute_phis(rp, rm, bp2i, bm2i, bpm2i, &
                                    bpi, bmi, bpmi, cosphim2, cosphim1, phim1, phim2)
                                call compute_phis(rm, rp, bm2i, bp2i, bpm2i, &
                                    bmi, bpi, bpmi, cosphip2, cosphip1, phip1, phip2)
                                lc(i) = 1.d0 - (Arc(c1, c2, phip1, phip2, rp, bpi) &
                                    + Arc(c1, c2, phim1, phim2, rm, bmi)) * of0
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
                                lc(i) = (2 * F(c1, c2, pi - thetam, 1.d0, 0.d0) &
                                  - Arc(c1, c2, -thetamstar, phim2, rm, bmi) &
                                  - Arc(c1, c2, phim1, thetamstar, rm, bmi) &
                                  - Arc(c1, c2, phip1, phip2, rp, bpi)) * of0
                                goto 1
                            end if
                        end if
                    else
                        if (bmi + rm .lt. 1.d0) then 
                            if (bpmi + rm .lt. rp) then
                                call compute_theta(rp, rp2, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                                call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costheta, theta)
                                lc(i) = 2 * (F(c1, c2, pi - theta, 1.d0, 0.d0) &
                                    - F(c1, c2, thetapstar, rp, bpi)) * of0
                                goto 1
                            else
                                call compute_phis(rp, rm, bp2i, bm2i, bpm2i,&
                                    bpi, bmi, bpmi, cosphim2, cosphim1, phim1, phim2)
                                call compute_phis(rm, rp, bm2i, bp2i, bpm2i, &
                                    bmi, bpi, bpmi, cosphip2, cosphip1, phip1, phip2)
                                call compute_theta(rp, rp2, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                                call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costhetap, thetap)
                                lc(i) = (2 * F(c1, c2, pi - thetap, 1.d0, 0.d0) &
                                  - Arc(c1, c2, -thetapstar, phip2, rp, bpi) &
                                  - Arc(c1, c2, phip1, thetapstar, rp, bpi) &
                                  - Arc(c1, c2, phim1, phim2, rm, bmi)) * of0
                                goto 1
                            end if
                        else
                            if (bpmi + rm .lt. rp) then
                                call compute_theta(rp, rp, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                                call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costheta, theta)
                                lc(i) = 2 * (F(c1, c2, pi - theta, 1.d0, 0.d0) &
                                    - F(c1, c2, thetapstar, rp, bpi)) * of0
                                goto 1
                            else
                                costhetapm = (bp2i + bm2i - bpm2i) / (2 * bpi * bmi)
                                cosphim = (bm2i + 1.d0 - rm**2.d0) / (2 * bmi)
                                cosphip = (bp2i + 1.d0 - rp**2.d0) / (2 * bpi)
                                if (abs(costhetapm - 1.d0) .lt. 1.d-5) then
                                    thetapm = 0.d0
                                else
                                    thetapm = Acos(costhetapm)
                                end if
                                phim = Acos(cosphim)
                                phip = Acos(cosphip)
                                if (thetapm + phim .lt. phip) then
                                    if ((bmi - rm) .lt. (bpi - rp)) then
                                        call compute_phis(rp, rm, bp2i, bm2i, bpm2i, &
                                            bpi, bmi, bpmi, cosphim2, cosphim1, phim1, phim2)
                                        call compute_phis(rm, rp, bm2i, bp2i, bpm2i, &
                                            bmi, bpi, bpmi, cosphip2, cosphip1, phip1, phip2)
                                        call compute_theta(rp, rp, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                                        call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costhetap, thetap)
                                        lc(i) = (2 * F(c1, c2, pi - thetap, 1.d0, 0.d0) &
                                          - Arc(c1, c2, -thetapstar, phip2, rp, bpi) &
                                          - Arc(c1, c2, phip1, thetapstar, rp, bpi) &
                                          - Arc(c1, c2, phim1, phim2, rm, bmi)) * of0
                                         goto 1
                                    else
                                        call compute_theta(rp, rp, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                                        call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costheta, theta)
                                        lc(i) = 2 * (F(c1, c2, pi - theta, 1.d0, 0.d0) &
                                          - F(c1, c2, thetapstar, rp, bpi)) * of0
                                        goto 1
                                    end if
                                else if (thetapm + phip .lt. phim) then
                                    if ((bpi - rp) .lt. (bmi - rm)) then
                                        call compute_phis(rp, rm, bp2i, bm2i, bpm2i, &
                                            bpi, bmi, bpmi, cosphim2, cosphim1, phim1, phim2)
                                        call compute_phis(rm, rp, bm2i, bp2i, bpm2i, &
                                            bmi, bpi, bpmi, cosphip2, cosphip1, phip1, phip2)
                                        call compute_theta(rm, rm, 1.d0, 1.d0, bmi, bm2i, costhetamstar, thetamstar)
                                        call compute_theta(1.d0, 1.d0, rm, rm2, bmi, bm2i, costhetam, thetam)
                                        lc(i) = (2 * F(c1, c2, pi - thetam, 1.d0, 0.d0) &
                                            - Arc(c1, c2, -thetamstar, phim2, rm, bmi) &
                                            - Arc(c1, c2, phim1, thetamstar, rm, bmi) & 
                                            - Arc(c1, c2, phip1, phip2, rp, bpi)) * of0
                                        goto 1
                                    else
                                        call compute_theta(rm, rm2, 1.d0, 1.d0, bmi, bm2i, costhetamstar, thetamstar)
                                        call compute_theta(1.d0, 1.d0, rm, rm2, bmi, bm2i, costheta, theta)
                                        lc(i) = 2 * (F(c1, c2, pi - theta, 1.d0, 0.d0) &
                                          - F(c1, c2, thetamstar, rm, bmi)) * of0
                                        goto 1
                                    end if
                                else
                                    costheta = (bpm2i + bm2i - bp2i) / (2 * bpmi * bmi)
                                    cosphim = (bpm2i + rm**2 - rp**2) / (2 * bpmi * rm)
                                    cosphi1 = Cos(Acos(costheta) - Acos(cosphim))
                                    cosphi2 = Cos(Acos(costheta) + Acos(cosphim))
                                    d1 = rm**2 + bm2i - 2 * rm * bmi * cosphi1
                                    d2 = rm**2 + bm2i - 2 * rm * bmi * cosphi2
                                    if (d1 .gt. 1.d0) then
                                        call compute_theta(rp, rp, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                                        call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costhetap, thetap)
                                        call compute_theta(rm, rm2, 1.d0, 1.d0, bmi, bm2i, costhetamstar, thetamstar)
                                        call compute_theta(1.d0, 1.d0, rm, rm2, bmi, bm2i, costhetam, thetam)
                                        theta = thetap + thetam
                                        lc(i) = 2 * (F(c1, c2, pi - theta, 1.d0, 0.d0) &
                                          - F(c1, c2, thetapstar, rp, bpi) & 
                                          - F(c1, c2, thetamstar, rm, bmi)) * of0
                                        goto 1
                                    else if (d2 .lt. 1.d0) then
                                        call compute_phis(rp, rm, bp2i, bm2i, bpm2i, &
                                            bpi, bmi, bpmi, cosphim2, cosphim1, phim1, phim2)
                                        call compute_phis(rm, rp, bm2i, bp2i, bpm2i, &
                                            bmi, bpi, bpmi, cosphip2, cosphip1, phip1, phip2)
                                        call compute_theta(rm, rm2, 1.d0, 1.d0, bmi, bm2i, costhetamstar, thetamstar)
                                        call compute_theta(1.d0, 1.d0, rm, rm2, bmi, bm2i, costhetam, thetam)
                                        call compute_theta(rp, rp, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                                        call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costhetap, thetap)
                                        theta = thetam + thetap
                                        lc(i) = (2 * F(c1, c2, pi - theta, 1.d0, 0.d0) &
                                            - Arc(c1, c2, -thetamstar, -phim1, rm, bmi) & 
                                            - Arc(c1, c2, -phim2, thetamstar, rm, bmi) &
                                            - Arc(c1, c2, phip1, thetapstar, rp, bpi) &
                                            - Arc(c1, c2, -thetapstar, phip2, rp, bpi)) * of0
                                        goto 1
                                else
                                        call compute_phis(rp, rm, bp2i, bm2i, bpm2i, &
                                            bpi, bmi, bpmi, cosphim2, cosphim1, phim1, phim2)
                                        call compute_phis(rm, rp, bm2i, bp2i, bpm2i, &
                                            bmi, bpi, bpmi, cosphip2, cosphip1, phip1, phip2)
                                        call compute_theta(rm, rm2, 1.d0, 1.d0, bmi, bm2i, costhetamstar, thetamstar)
                                        call compute_theta(1.d0, 1.d0, rm, rm2, bmi, bm2i, costhetam, thetam)
                                        call compute_theta(rp, rp, 1.d0, 1.d0, bpi, bp2i, costhetapstar, thetapstar)
                                        call compute_theta(1.d0, 1.d0, rp, rp2, bpi, bp2i, costhetap, thetap)
                                        theta = 0.5 * (thetap + thetam + thetapm)
                                        lc(i) = (2 * F(c1,c2, pi - theta, 1.d0, 0.d0) & 
                                            - Arc(c1, c2, -thetamstar, -phim1, rm, bmi) & 
                                            - Arc(c1, c2, phip1, thetapstar, rp, bpi)) * of0
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

real*8 function Arc(c1, c2, phi1, phi2, r, b)

    real*8 :: phi1, phi2, r, b, c1, c2
        
    if (phi1 < 0) then
        if (phi2 < 0) then
            if (phi2 < phi1) then
                Arc = 2 * F(c1, c2, pi, r, b) + F(c1, c2, -phi1, r, b) - F(c1, c2, -phi2, r, b)
                return
            else
                Arc = -F(c1, c2, -phi2, r, b) + F(c1, c2, -phi1, r, b)
                return
            end if
        else
            Arc = F(c1, c2, phi2, r, b) + F(c1, c2, -phi1, r, b)
            return
        end if
    else
        if (phi2 < 0) then
            Arc = 2 * F(c1, c2, pi, r, b) - F(c1, c2, phi1, r, b) - F(c1, c2, -phi2, r, b)
            return
        else
            if (phi2 < phi1) then
                Arc = 2 * F(c1, c2, pi, r, b) + F(c1, c2, phi2, r, b) - F(c1, c2, phi1, r, b)
                return
            else
                Arc = F(c1, c2, phi2, r, b) - F(c1, c2, phi1, r, b)
                return
            end if
        end if
    end if
    
    return
    
end function

real*8 function F(c1, c2, phi, r, b)

    real*8 :: c1, c2, phi, r, b, Fc, Fq, Fl
    real*8 :: o, pc, ac, bc
    real*8 :: alpha, beta, gamma, d, s, n, m, x, y, ome, tans, sphi
    real*8 :: s2phi, s3phi, cphi, s4phi, r2, b2, bmr2, tanphi2, sinphi2
    real*8 :: ellipf, ellipe, ellippi, ellipf_tmp, ome_tmp, tans_tmp
    
    sphi = Sin(phi)
    cphi = Cos(phi)
    s2phi = 2 * sphi * cphi
    s3phi = 3 * sphi - 4 * sphi * sphi * sphi
    s4phi = 4 * sphi * cphi * cphi * cphi - 4 * sphi * sphi * sphi * cphi
    r2 = r * r
    b2 = b * b
    bmr2 = (b - r)**2.d0
    tanphi2 = sphi / (cphi + 1.d0)
    sinphi2 = Sin(phi * 0.5)
    
    Fc = 0.5 * r * (r * phi - b * sphi)
    
    Fq = -(r * (4 * b * (2 * b2 + 9 * r2) * sphi &
             - 4 * r * (3 * (2 * b2 + r2) * phi + b * r * s3phi) &
             + r * r2 * s4phi)) / 48.d0
    
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
            ellipf = el1(Tan(s), Sqrt(1.d0 - m))
            o = 1.d0
            ellipe = el2(Tan(s), Sqrt(1.d0 - m), o, 1.d0 - m)
            Fl = alpha * ellipe + beta * ellipf + d
            goto 2
        end if
    else if (b + r == 1.d0) then
        y = Sqrt(b * r)
        Fl = phi * o3 * 0.5 + Atan((2 * r - 1.d0) / tanphi2) * o3 &
                + Atan(2 * y * sinphi2 / (1.d0 - 2 * r)) * o3 &
                + pi *  Sqrt(1.d0 - 4 * r * b) / (2 * r - 1.d0) * o3 * 0.5 &
                + 2.d0 * o9 * y * (2 * r * (5 * r - 2.d0) &
                - 2 * b * r * cphi - 3.d0) * sinphi2
        goto 2
    else if (b + r .gt. 1.d0) then
        
        x = Sqrt(1.d0 - bmr2)
        y = Sqrt(b * r)
        alpha = 2 * y * (7 * r2 + b2 - 4.d0) * o9
        beta = 3.d0 + 2*r*(b2 * b + 5 * b2 * r + 3*r*(-2.d0 + r2) + b*(-4.d0 + 7*r2))
        beta = -beta * o9 * 0.5 / y
        gamma = (b + r) * o3 / (2 * (b - r) * y)
        m = (1.d0 - bmr2) / (4 * r * b)
        ome = Sqrt(1.d0 - m)
        if (abs(sinphi2 / Sqrt(m) - 1.D0) .lt. 1.D-7) then
            d = phi * o3 * 0.5 - Atan((b + r) / (b - r) * tanphi2) * o3
            n = 1.d0 - 1.d0 / bmr2
            pc = 1.d0
            ac = 1.d0
            bc = 1.d0
            o = 1.d0
            ome_tmp = ome
            ellipf = cel(ome_tmp, pc, ac, bc)
            pc = 1.d0
            ac = 1.d0
            ome_tmp = ome
            ellipe = cel(ome_tmp, pc, ac, 1.d0 - m)
            ac = 1.d0
            bc = 1.d0
            ome_tmp = ome
            ellippi = cel(ome_tmp, 1.d0 - n, ac, bc)
        else
            d = phi * o3 * 0.5 - Atan((b + r) / (b - r) * tanphi2) * o3 &
                - (2 * b * r * o9) * sphi &
                * Sqrt(1.d0 - b2 - r2 + 2 * b * r * cphi)
            s = Asin(sinphi2 / Sqrt(m))
            tans = Tan(s)
            n = 1.d0 - 1.d0 / bmr2
            tans_tmp = tans
            ome_tmp = ome
            ellipf = el1(tans_tmp, ome_tmp)
            o = 1.d0
            tans_tmp = tans
            ome_tmp = ome
            ellipe = el2(tans_tmp, ome_tmp, o, 1.d0 - m)
            tans_tmp = tans
            ome_tmp = ome
            ellippi = el3(tans_tmp, ome_tmp, 1.d0 - n)
        end if
        
        Fl = alpha * ellipe + beta * ellipf + gamma * ellippi + d
        goto 2
    else
        x = Sqrt(1.d0 - bmr2)
        alpha = (7 * r2 + b2 - 4.d0) * x * o9
        beta = (r**4.d0 + b**4.d0 + r2 - b2 * (5.d0 + 2 * r2) + 1.d0) / (9.d0 * x)
        gamma = (b + r) / (b - r) / (3.d0 * x)
        n = - 4 * r * b / bmr2
        m = 4 * r * b / (1.d0 - bmr2)
        ome = Sqrt(1.d0 - m)
        if (abs(abs(phi) - pi) .lt. 1.D-5) then
        
            d = pi * 0.5 * o3 * (1.d0 - Sign(1.d0, b - r))
        
            pc = 1.d0
            ac = 1.d0
            bc = 1.d0
            o = 1.d0
            ome_tmp = ome
            ellipf = cel(ome_tmp, pc, ac, bc)
            pc = 1.d0
            ac = 1.d0
            ome_tmp = ome
            ellipe = cel(ome_tmp, pc, ac, 1.d0 - m)
            ac = 1.d0
            bc = 1.d0
            ome_tmp = ome
            ellippi = cel(ome_tmp, 1.d0 - n, ac, bc) 
        else
            d = phi * o3 * 0.5 - Atan((b + r) / (b - r) * tanphi2) * o3 &
                - (2 * b * r * o9) * sphi &
                * Sqrt(1.d0 - b2 - r2 + 2 * b * r * cphi)
            ome_tmp = ome
            tans_tmp = tanphi2
            ellipf = el1(tans_tmp, ome_tmp)
            o = 1.d0
            ome_tmp = ome
            tans_tmp = tanphi2
            ellipe = el2(tans_tmp, ome_tmp, o, 1.d0 - m)
            ome_tmp = ome
            tans_tmp = tanphi2
            ellippi = el3(tans_tmp, ome_tmp, 1.d0 - n)
        end if
        Fl = alpha * ellipe + beta * ellipf + gamma * ellippi + d
        goto 2
    end if
    
2   F = (1.d0 - c1 - 2 * c2) * Fc + (c1 + 2 * c2) * Fl + c2 * Fq
    
    return
end function

real*8 function F_const(phi, r, b)

    real*8 :: phi, r, b
    
    F_const = 0.5 * r * (r * phi - b * Sin(phi))
    return
    
end function

real*8 function F_lin(phi, r, b)

    real*8 :: phi, r, b
    real*8 :: o
    real*8 :: alpha, beta, gamma, d, s, n, m, x
    real*8 :: ellipf, ellipe, ellippi, ellipf_tmp
    
    if (b == 0.d0) then
        if (r == 1.d0) then
            F_lin = phi * o3
            return
        else
            F_lin = phi * (1.d0 - (1.d0 - r * r) ** (1.5)) * o3
            return
        end if
    else if (b == r) then
        if (r == 0.5) then
            s = phi * 0.5
            F_lin = phi * o3 * 0.5 + o3 * Sin(s) &
                   * (1.d0 - Sin(s)**2.d0 * o3)
            return
        else
            alpha = 4 * (2 * r * r - 1.d0) * o9
            beta = (1.d0 - 4 * r * r) * o9
            d = phi * o3 * 0.5 - 2 * r * r * Sin(phi) & 
                * Sqrt(1.d0 + 2 * r * r * (Cos(phi)-1.d0)) * o9
            s = phi * 0.5
            m = 4 * r * r
            ellipf = el1(Tan(s), Sqrt(1.d0 - m))
            o = 1.d0
            ellipe = el2(Tan(s), Sqrt(1.d0 - m), o, 1.d0 - m)
            F_lin = alpha * ellipe + beta * ellipf + d
            return
        end if
    else if (b + r == 1.d0) then
        s = phi * 0.5
        F_lin = phi * o3 * 0.5 + Atan((2 * r - 1.d0) / Tan(s)) * o3 &
                + Atan(2 * Sqrt(r * b) * Sin(s) / (1.d0 - 2 * r)) * o3 &
                + pi *  Sqrt(1.d0 - 4 * r * b) / (2 * r - 1.d0) * o3 * 0.5 &
                + 2.d0 * o9 * Sqrt(b * r) * (2 * r * (5 * r - 2.d0) &
                - 2 * b * r * Cos(phi) - 3.d0) * Sin(s)
        return
    else if (b + r .gt. 1.d0) then
        
        x = Sqrt(1.d0 - (b - r)**2.d0)
        s = phi * 0.5
        alpha = 2 * Sqrt(b * r) * (7 * r * r + b * b - 4.d0) * o9
        beta = 3.d0 + 2*r*(b**3.d0 + 5*b*b*r + 3*r*(-2.d0 + r*r) + b*(-4.d0 + 7*r*r))
        beta = -beta * o9 * 0.5 / Sqrt(b * r)
        gamma = (b + r) * o3 / (2 * (b - r) * Sqrt(b * r))
        d = phi * o3 * 0.5 - Atan((b + r) / (b - r) * Tan(s)) * o3 &
                - (2 * b * r * o9) * Sin(phi) &
                * Sqrt(1.d0 - b * b - r * r + 2 * b * r * Cos(phi))
        
        m = (1.d0 - (r - b)**2.d0) / (4 * r * b)
        if (abs(Sin(s) / Sqrt(m) - 1.D0) .lt. 1.D-7) then
            d = phi * o3 * 0.5 - Atan((b + r) / (b - r) * Tan(s)) * o3
            s = pi * Sign(1.d0, s) * 0.5
        else
            s = Asin(Sin(s) / Sqrt(m))
        end if

        n = 1.d0 - 1.d0 / (b - r)**2.d0
        ellipf = el1(Tan(s), Sqrt(1.d0 - m))
        o = 1.d0
        ellipe = el2(Tan(s), Sqrt(1.d0 - m), o, 1.d0 - m)
        ellippi = el3(Tan(s), Sqrt(1.d0 - m), 1.d0 - n)
        
        F_lin = alpha * ellipe + beta * ellipf + gamma * ellippi + d
        return
    else
        x = Sqrt(1.d0 - (b - r)**2.d0)
        s = phi * 0.5
        alpha = (7 * r * r + b * b - 4.d0) * x * o9
        beta = (r**4.d0 + b**4.d0 + r * r - b * b * (5.d0 + 2 * r * r) + 1.d0) / (9.d0 * x)
        gamma = (b + r) / (b - r) / (3.d0 * x)
        d = phi * o3 * 0.5 - Atan((b + r) / (b - r) * Tan(s)) * o3 &
                - (2 * b * r * o9) * Sin(phi) &
                * Sqrt(1.d0 - b * b - r * r + 2 * b * r * Cos(phi))
        n = - 4 * r * b / (b - r)**2.d0
        m = 4 * r * b / (1.d0 - (r - b)**2.d0)
        ellipf = el1(Tan(s), Sqrt(1.d0 - m))
        o = 1.d0
        ellipe = el2(Tan(s), Sqrt(1.d0 - m), o, 1.d0 - m)
        ellippi = el3(Tan(s), Sqrt(1.d0 - m), 1.d0 - n)
        F_lin = alpha * ellipe + beta * ellipf + gamma * ellippi + d
        return
    end if
        
    F_lin = 0.d0
    return
end function

real*8 function F_quad(phi, r, b)

    real*8 :: phi, r, b
    
    F_quad = -(r*(4*b*(2*b**2.d0 + 9*r**2.d0)*Sin(phi) &
             - 4*r*(3*(2*b**2.d0 + r**2.d0)*phi + b*r*Sin(3*phi)) &
             + r**3.d0*Sin(4*phi))) / 48.d0
    
    return
end function

subroutine F_const_wrapper(phi, r, b, res) bind(C, name="F_const")

    real*8, bind(C) :: phi, r, b
    real*8, bind(C), intent(out) :: res
    res = F_const(phi, r, b)
    
end

subroutine F_lin_wrapper(phi, r, b, res) bind(C, name="F_lin")

    real*8, bind(C) :: phi, r, b
    real*8, bind(C), intent(out) :: res
    res = F_lin(phi, r, b)
    
end

subroutine F_quad_wrapper(phi, r, b, res) bind(C, name="F_quad")

    real*8, bind(C) :: phi, r, b
    real*8, bind(C), intent(out) :: res
    res = F_quad(phi, r, b)
    
end

subroutine F_wrapper(c1, c2, phi, r, b, res) bind(C, name="F")

    real*8, bind(C) :: c1, c2, phi, r, b
    real*8, bind(C), intent(out) :: res
    res = F(c1, c2, phi, r, b)

end

subroutine Arc_wrapper(c1, c2, phi1, phi2, r, b, res) bind(C, name="Arc")

    real*8, bind(C) :: phi1, phi2, r, b, c1, c2
    real*8, bind(C), intent(out) :: res
    res = Arc(c1, c2, -pi*0.999999, pi*0.999999, r, b)
    
end

subroutine ellipf(phi, m, e) bind(C, name="f_burl")

    real*8, bind(C) :: phi, m
    real*8, bind(C), intent(out) :: e
    e = el1(Tan(phi), Sqrt(1.d0 - m))
end
    
subroutine ellipe(phi, m, e) bind(C, name="e_burl")

    real*8, bind(C) :: phi, m
    real*8, bind(C), intent(out) :: e
    real*8 :: o
    o = 1.d0
    e = el2(Tan(phi), Sqrt(1.d0 - m), o, 1.d0 - m)
end
    
subroutine ellipp(phi, n, m, e) bind(C, name="p_burl")

    real*8, bind(C) :: phi, m, n
    real*8, bind(C), intent(out) :: e
    e = el3(Tan(phi), Sqrt(1.d0 - m), 1.d0 - n)
end

subroutine ellipkc(kc, e) bind(C, name="kc_burl")
   
   real*8, bind(C) :: kc
   real*8, bind(C), intent(out) :: e
   real*8 :: p, a, b
   p = 1.d0
   a = 1.d0
   b = 1.d0
   e = cel(kc, p, a, b)
end

subroutine ellipec(kc, e) bind(C, name="ec_burl")
   
   real*8, bind(C) :: kc
   real*8, bind(C), intent(out) :: e
   real*8 :: p, a, b
   p = 1.d0
   a = 1.d0
   b = kc * kc
   e = cel(kc, p, a, b)
end

subroutine ellippc(kc, p, e) bind(C, name="pc_burl")
   
   real*8, bind(C) :: kc, p
   real*8, bind(C), intent(out) :: e
   real*8 :: a, b
   a = 1.d0
   b = 1.d0
   e = cel(kc, p, a, b)
end

end module phot