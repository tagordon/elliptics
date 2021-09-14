module phot
use iso_c_binding
use ellip

implicit none

real*8, parameter :: pi = 3.14159265358979323846, pilims = 3.1415926
real*8, parameter :: o3 = 0.33333333333333333333, o9 = 0.1111111111111111111

contains

subroutine compute_phis(rp, rm, bp2, bm2, bpm2, cosphi2, cosphi1, phi1, phi2)

    real*8 :: rp, rm, bp2, bm2, bpm2
    real*8, intent(out) :: cosphi1, cosphi2, phi1, phi2
    real*8 :: delta, gamma, s, t, tmp
    
    real*8 :: bp, bm, bpm
    bp = Sqrt(bp2)
    bm = Sqrt(bm2)
    bpm = Sqrt(bpm2)
    
    delta = bm * bpm2 * rm &
          * Sqrt(((bm - bp - bpm) * (bm + bp - bpm) * (bm - bp + bpm) &
          * (bm + bp + bpm) * (bpm - rm - rp) * (bpm + rm - rp) &
          * (bpm - rm + rp) * (bpm + rm + rp)) / (bm2 * bpm2**2 * rm**2))
          
    gamma = bm2 * (bpm2 + rm**2 - rp**2) + (-bp2 + bpm2)*(bpm2 + rm**2 - rp**2)
    cosphi2 = (gamma - delta) / (4. * bm * bpm2 * rm)
    cosphi1 = (gamma + delta) / (4. * bm * bpm2 * rm)
    phi1 = Acos(cosphi1)
    phi2 = Acos(cosphi2)
    
    s = (bm2**2 + (bp2 - bpm2)**2. & 
      - (bm2*(bpm2**2. + 2. * bp2 * rm**2 & 
      + rm**4 - 2. * (bpm2 + rm**2.) * rp**2. + rp**4)) / rm**2) & 
      / (4. * bm2 * bpm**2)
    s = Sign(1.d0, s)
    t = bpm**2.d0 - (-(bm2 * rm) + bp2 * rm - bm * rm**2. + bm * rp**2.d0) & 
      / (bm + rm)
    t = Sign(1.d0, t)
    phi1 = s * t * phi1
    phi2 = - t * phi2
end

! angle is with respect to the r1 circle's center
subroutine compute_theta(r1, r2, b, costheta, theta)

    real*8 :: r1, r2, b, costheta, theta
    
    costheta = (r1**2.d0 + b**2.d0 - r2**2.d0) / (2.d0 * r1 * b)
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
    real*8 :: d1, d2, a
    real*8 :: f0, of0
    
    real*8 :: bp(j), bm(j), bpm(j)
    bp = Sqrt(bp2)
    bm = Sqrt(bm2)
    bpm = Sqrt(bpm2)
    
    f0 = (1.d0 - c1 - 2 * c2) * pi + (c1 + 2 * c2) * 2 * pi * o3 + c2 * 0.5 * pi
    of0 = 1.d0 / f0
    lc = 1.d0
    
    do i=1,j,1
        
        if ((bp(i) .gt. rp + 1.d0) .AND. (bm(i) .gt. rm + 1.d0)) then
            lc(i) = 1.d0
            goto 1
        
        else if (bpm(i) .gt. rp + rm) then
            if (bp(i) .gt. rp + 1.d0) then
                if (bm(i) .gt. rm + 1.d0) then
                    lc(i) = 1.d0
                    goto 1
                else
                    if (bm(i) + rm .lt. 1.d0) then
                        lc(i) = 1.d0 - Arc(c1, c2, -pilims, pilims, rm, bm(i)) * of0
                        goto 1
                    else
                        call compute_theta(rm, 1.d0, bm(i), costhetamstar, thetamstar)
                        call compute_theta(1.d0, rm, bm(i), costheta, theta)
                        lc(i) = (Arc(c1, c2, theta - pilims, pilims - theta, 1.d0, 0.d0) &
                              - Arc(c1, c2, -thetamstar, thetamstar, rm, bm(i))) * of0
                        goto 1
                    end if
                end if
            else
                if (bm(i) .gt. rm + 1.d0) then
                    if (bp(i) + rp .lt. 1.d0) then
                        lc(i) = 1.d0 - Arc(c1, c2, -pilims, pilims, rp, bp(i)) * of0
                        goto 1
                    else
                        call compute_theta(rp, 1.d0, bp(i), costhetapstar, thetapstar)
                        call compute_theta(1.d0, rp, bp(i), costheta, theta)
                        lc(i) = (Arc(c1, c2, theta - pilims, pilims - theta, 1.d0, 0.d0) &
                              - Arc(c1, c2, -thetapstar, thetapstar, rp, bp(i))) * of0
                        goto 1
                    end if
                else
                    if (bp(i) + rp .lt. 1.d0) then
                        if (bm(i) + rm .lt. 1.d0) then
                            lc(i) = 1.d0 - (Arc(c1, c2, -pilims, pilims, rm, bm(i)) &
                                  + Arc(c1, c2, -pilims, pilims, rp, bp(i))) * of0
                            goto 1
                        else
                            call compute_theta(rm, 1.d0, bm(i), costhetamstar, thetamstar)
                            call compute_theta(1.d0, rm, bm(i), costheta, theta)
                            lc(i) = (Arc(c1, c2, theta - pilims, pilims - theta, 1.d0, 0.d0) &
                                  - Arc(c1, c2, -thetamstar, thetamstar, rm, bm(i)) &
                                  - Arc(c1, c2, -pilims, pilims, rp, bp(i))) * of0
                            goto 1
                        end if
                    else
                        if (bm(i) + rm .lt. 1.d0) then
                            call compute_theta(rp, 1.d0, bp(i), costhetapstar, thetapstar)
                            call compute_theta(1.d0, rp, bp(i), costheta, theta)
                            lc(i) = (Arc(c1, c2, theta - pilims, pilims - theta, 1.d0, 0.d0) &
                                  - Arc(c1, c2, -thetapstar, thetapstar, rp, bp(i)) & 
                                  - Arc(c1, c2, -pilims, pilims, rm, bm(i))) * of0
                            goto 1
                        else
                            call compute_theta(rp, 1.d0, bp(i), costhetapstar, thetapstar)
                            call compute_theta(1.d0, rp, bp(i), costhetap, thetap)
                            call compute_theta(rm, 1.d0, bm(i), costhetamstar, thetamstar)
                            call compute_theta(1.d0, rm, bm(i), costhetam, thetam)
                            theta = thetap + thetam
                            lc(i) = (Arc(c1, c2, theta - pilims, pilims - theta, 1.d0, 0.d0) &
                                  - Arc(c1, c2, -thetapstar, thetapstar, rp, bp(i)) & 
                                  - Arc(c1, c2, -thetamstar, thetamstar, rm, bm(i))) * of0
                            goto 1
                        end if
                    end if
                end if
            end if
        else
            if (bp(i) .gt. rp + 1.d0) then
                if (bm(i) .gt. rm + 1.d0) then
                    lc(i) = 1.d0
                    goto 1
                else
                    if (bm(i) + rm .lt. 1.d0) then
                        lc(i) = 1.d0 - Arc(c1, c2, -pilims, pilims, rm, bm(i)) * of0
                        goto 1
                    else
                        call compute_theta(rm, 1.d0, bm(i), costhetamstar, thetamstar)
                        call compute_theta(1.d0, rm, bm(i), costheta, theta)
                        lc(i) = (Arc(c1, c2, theta - pilims, pilims - theta, 1.d0, 0.d0) &
                              - Arc(c1, c2, -thetamstar, thetamstar, rm, bm(i))) * of0
                        goto 1
                    end if
                end if
            else
                if (bm(i) .gt. rm + 1.d0) then
                    if (bp(i) + rp .lt. 1.d0) then
                        lc(i) = 1.d0 - Arc(c1, c2, -pilims, pilims, rp, bp(i)) * of0
                        goto 1
                    else
                        call compute_theta(rp, 1.d0, bp(i), costhetapstar, thetapstar)
                        call compute_theta(1.d0, rp, bp(i), costheta, theta)
                        lc(i) = (Arc(c1, c2, theta - pilims, pilims - theta, 1.d0, 0.d0) &
                              - Arc(c1, c2, -thetapstar, thetapstar, rp, bp(i))) * of0
                        goto 1
                    end if
                else
                    if (bp(i) + rp .lt. 1.d0) then
                        if (bm(i) + rm .lt. 1.d0) then
                            if (bpm(i) + rm .lt. rp) then
                                lc(i) = 1.d0 - Arc(c1, c2, -pilims, pilims, rp, bp(i)) * of0
                                goto 1
                            else
                                call compute_phis(rp, rm, bp2(i), bm2(i), bpm2(i), cosphim2, cosphim1, phim1, phim2)
                                call compute_phis(rm, rp, bm2(i), bp2(i), bpm2(i), cosphip2, cosphip1, phip1, phip2)
                                lc(i) = 1.d0 - (Arc(c1, c2, phip1, phip2, rp, bp(i)) &
                                      + Arc(c1, c2, phim1, phim2, rm, bm(i))) * of0
                                goto 1
                            end if
                        else
                            if (bpm(i) + rm .lt. rp) then
                                goto 1
                            else
                                call compute_phis(rp, rm, bp2(i), bm2(i), bpm2(i), cosphim2, cosphim1, phim1, phim2)
                                call compute_phis(rm, rp, bm2(i), bp2(i), bpm2(i), cosphip2, cosphip1, phip1, phip2)
                                call compute_theta(rm, 1.d0, bm(i), costhetamstar, thetamstar)
                                call compute_theta(1.d0, rm, bm(i), costhetam, thetam)
                                lc(i) = (Arc(c1, c2, thetam - pilims, pilims - thetam, 1.d0, 0.d0) &
                                      - Arc(c1, c2, -thetamstar, phim2, rm, bm(i)) &
                                      - Arc(c1, c2, phim1, thetamstar, rm, bm(i)) &
                                      - Arc(c1, c2, phip1, phip2, rp, bp(i))) * of0
                                goto 1
                            end if
                        end if
                    end if
                                
            
                    if (bp(i) + rp .lt. 1.d0) then
                        if (bm(i) + rm .lt. 1.d0) then
                            if (bpm(i) + rm .lt. rp) then
                                lc(i) = 1.d0 - Arc(c1, c2, -pilims, pilims, rp, bp(i)) * of0
                                goto 1
                            else
                                call compute_phis(rp, rm, bp2(i), bm2(i), bpm2(i), cosphim2, cosphim1, phim1, phim2)
                                call compute_phis(rm, rp, bm2(i), bp2(i), bpm2(i), cosphip2, cosphip1, phip1, phip2)
                                lc(i) = 1.d0 - (Arc(c1, c2, phip1, phip2, rp, bp(i)) + Arc(c1, c2, phim1, phim2, rm, bm(i))) * of0
                                goto 1
                            end if
                        else 
                            if (bpm(i) + rm .lt. rp) then
                                goto 1
                            else
                                call compute_phis(rp, rm, bp2(i), bm2(i), bpm2(i), cosphim2, cosphim1, phim1, phim2)
                                call compute_phis(rm, rp, bm2(i), bp2(i), bpm2(i), cosphip2, cosphip1, phip1, phip2)
                                call compute_theta(rm, 1.d0, bm(i), costhetamstar, thetamstar)
                                call compute_theta(1.d0, rm, bm(i), costhetam, thetam)
                                lc(i) = (Arc(c1, c2, thetam - pilims, pilims - thetam, 1.d0, 0.d0) &
                                  - Arc(c1, c2, -thetamstar, phim2, rm, bm(i)) &
                                  - Arc(c1, c2, phim1, thetamstar, rm, bm(i)) &
                                  - Arc(c1, c2, phip1, phip2, rp, bp(i))) * of0
                                goto 1
                            end if
                        end if
                    else
                        if (bm(i) + rm .lt. 1.d0) then 
                            if (bpm(i) + rm .lt. rp) then
                                call compute_theta(rp, 1.d0, bp(i), costhetapstar, thetapstar)
                                call compute_theta(1.d0, rp, bp(i), costheta, theta)
                                lc(i) = (Arc(c1, c2, theta - pilims, pilims - theta, 1.d0, 0.d0) &
                                    - Arc(c1, c2, -thetapstar, thetapstar, rp, bp(i))) * of0
                                goto 1
                            else
                                call compute_phis(rp, rm, bp2(i), bm2(i), bpm2(i), cosphim2, cosphim1, phim1, phim2)
                                call compute_phis(rm, rp, bm2(i), bp2(i), bpm2(i), cosphip2, cosphip1, phip1, phip2)
                                call compute_theta(rp, 1.d0, bp(i), costhetapstar, thetapstar)
                                call compute_theta(1.d0, rp, bp(i), costhetap, thetap)
                                lc(i) = (Arc(c1, c2, thetap - pilims, pilims - thetap, 1.d0, 0.d0) &
                                  - Arc(c1, c2, -thetapstar, phip2, rp, bp(i)) &
                                  - Arc(c1, c2, phip1, thetapstar, rp, bp(i)) &
                                  - Arc(c1, c2, phim1, phim2, rm, bm(i))) * of0
                                goto 1
                            end if
                        else
                            if (bpm(i) + rm .lt. rp) then
                                call compute_theta(rp, 1.d0, bp(i), costhetapstar, thetapstar)
                                call compute_theta(1.d0, rp, bp(i), costheta, theta)
                                lc(i) = (Arc(c1, c2, theta - pilims, pilims - theta, 1.d0, 0.d0) &
                                    - Arc(c1, c2, -thetapstar, thetapstar, rp, bp(i))) * of0
                                goto 1
                            else
                                costhetapm = (bp2(i) + bm2(i) - bpm2(i)) / (2 * bp(i) * bm(i))
                                cosphim = (bm2(i) + 1.d0 - rm**2.d0) / (2 * bm(i))
                                cosphip = (bp2(i) + 1.d0 - rp**2.d0) / (2 * bp(i))
                                if (abs(costhetapm - 1.d0) .lt. 1.d-5) then
                                    thetapm = 0.d0
                                else
                                    thetapm = Acos(costhetapm)
                                end if
                                phim = Acos(cosphim)
                                phip = Acos(cosphip)
                                if (thetapm + phim .lt. phip) then
                                    if ((bm(i) - rm) .lt. (bp(i) - rp)) then
                                        call compute_phis(rp, rm, bp2(i), bm2(i), bpm2(i), cosphim2, cosphim1, phim1, phim2)
                                        call compute_phis(rm, rp, bm2(i), bp2(i), bpm2(i), cosphip2, cosphip1, phip1, phip2)
                                        call compute_theta(rp, 1.d0, bp(i), costhetapstar, thetapstar)
                                        call compute_theta(1.d0, rp, bp(i), costhetap, thetap)
                                        lc(i) = (Arc(c1, c2, thetap - pilims, pilims - thetap, 1.d0, 0.d0) &
                                          - Arc(c1, c2, -thetapstar, phip2, rp, bp(i)) &
                                          - Arc(c1, c2, phip1, thetapstar, rp, bp(i)) &
                                          - Arc(c1, c2, phim1, phim2, rm, bm(i))) * of0
                                         goto 1
                                    else
                                        call compute_theta(rp, 1.d0, bp(i), costhetapstar, thetapstar)
                                        call compute_theta(1.d0, rp, bp(i), costheta, theta)
                                        lc(i) = (Arc(c1, c2, theta - pilims, pilims - theta, 1.d0, 0.d0) &
                                          - Arc(c1, c2, -thetapstar, thetapstar, rp, bp(i))) * of0
                                        goto 1
                                    end if
                                else if (thetapm + phip .lt. phim) then
                                    if ((bp(i) - rp) .lt. (bm(i) - rm)) then
                                        call compute_phis(rp, rm, bp2(i), bm2(i), bpm2(i), cosphim2, cosphim1, phim1, phim2)
                                        call compute_phis(rm, rp, bm2(i), bp2(i), bpm2(i), cosphip2, cosphip1, phip1, phip2)
                                        call compute_theta(rm, 1.d0, bm(i), costhetamstar, thetamstar)
                                        call compute_theta(1.d0, rm, bm(i), costhetam, thetam)
                                        lc(i) = (Arc(c1, c2, thetam - pilims, pilims - thetam, 1.d0, 0.d0) &
                                            - Arc(c1, c2, -thetamstar, phim2, rm, bm(i)) &
                                            - Arc(c1, c2, phim1, thetamstar, rm, bm(i)) & 
                                            - Arc(c1, c2, phip1, phip2, rp, bp(i))) * of0
                                        goto 1
                                    else
                                        call compute_theta(rm, 1.d0, bm(i), costhetamstar, thetamstar)
                                        call compute_theta(1.d0, rm, bm(i), costheta, theta)
                                        lc(i) = (Arc(c1, c2, theta - pilims, pilims - theta, 1.d0, 0.d0) &
                                          - Arc(c1, c2, -thetamstar, thetamstar, rm, bm(i))) * of0
                                        goto 1
                                    end if
                                else
                                    costheta = (bpm2(i) + bm2(i) - bp2(i)) / (2 * bpm(i) * bm(i))
                                    cosphim = (bpm2(i) + rm**2 - rp**2) / (2 * bpm(i) * rm)
                                    cosphi1 = Cos(Acos(costheta) - Acos(cosphim))
                                    cosphi2 = Cos(Acos(costheta) + Acos(cosphim))
                                    d1 = rm**2 + bm2(i) - 2 * rm * bm(i) * cosphi1
                                    d2 = rm**2 + bm2(i) - 2 * rm * bm(i) * cosphi2
                                    if (d1 .gt. 1.d0) then
                                        call compute_theta(rp, 1.d0, bp(i), costhetapstar, thetapstar)
                                        call compute_theta(1.d0, rp, bp(i), costhetap, thetap)
                                        call compute_theta(rm, 1.d0, bm(i), costhetamstar, thetamstar)
                                        call compute_theta(1.d0, rm, bm(i), costhetam, thetam)
                                        theta = thetap + thetam
                                        lc(i) = (Arc(c1, c2, theta - pilims, pilims - theta, 1.d0, 0.d0) &
                                          - Arc(c1, c2, -thetapstar, thetapstar, rp, bp(i)) & 
                                          - Arc(c1, c2, -thetamstar, thetamstar, rm, bm(i))) * of0
                                        goto 1
                                    else if (d2 .lt. 1.d0) then
                                        call compute_phis(rp, rm, bp2(i), bm2(i), bpm2(i), cosphim2, cosphim1, phim1, phim2)
                                        call compute_phis(rm, rp, bm2(i), bp2(i), bpm2(i), cosphip2, cosphip1, phip1, phip2)
                                        call compute_theta(rm, 1.d0, bm(i), costhetamstar, thetamstar)
                                        call compute_theta(1.d0, rm, bm(i), costhetam, thetam)
                                        call compute_theta(rp, 1.d0, bp(i), costhetapstar, thetapstar)
                                        call compute_theta(1.d0, rp, bp(i), costhetap, thetap)
                                        theta = thetam + thetap
                                        lc(i) = (Arc(c1, c2, theta - pilims, pilims - theta, 1.d0, 0.d0) &
                                            - Arc(c1, c2, -thetamstar, -phim1, rm, bm(i)) & 
                                            - Arc(c1, c2, -phim2, thetamstar, rm, bm(i)) &
                                            - Arc(c1, c2, phip1, thetapstar, rp, bp(i)) &
                                            - Arc(c1, c2, -thetapstar, phip2, rp, bp(i))) * of0
                                        goto 1
                                else
                                        call compute_phis(rp, rm, bp2(i), bm2(i), bpm2(i), cosphim2, cosphim1, phim1, phim2)
                                        call compute_phis(rm, rp, bm2(i), bp2(i), bpm2(i), cosphip2, cosphip1, phip1, phip2)
                                        call compute_theta(rm, 1.d0, bm(i), costhetamstar, thetamstar)
                                        call compute_theta(1.d0, rm, bm(i), costhetam, thetam)
                                        call compute_theta(rp, 1.d0, bp(i), costhetapstar, thetapstar)
                                        call compute_theta(1.d0, rp, bp(i), costhetap, thetap)
                                        theta = 0.5 * (thetap + thetam + thetapm)
                                        lc(i) = (Arc(c1,c2, theta - pilims, pilims - theta, 1.d0, 0.d0) & 
                                            - Arc(c1, c2, -thetamstar, -phim1, rm, bm(i)) & 
                                            - Arc(c1, c2, phip1, thetapstar, rp, bp(i))) * of0
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
                Arc = 2 * F(c1, c2, pilims, r, b) + F(c1, c2, -phi1, r, b) - F(c1, c2, -phi2, r, b)
                return
            else
                Arc = -F(c1, c2, -phi2, r, b) + F(c1, c2, -phi1, r, b)
                return
            end if
        else
            if (phi1 == -phi2) then
                Arc = 2 * F(c1, c2, phi2, r, b)
                return
            else
                Arc = F(c1, c2, phi2, r, b) + F(c1, c2, -phi1, r, b)
                return
            end if
        end if
    else
        if (phi2 < 0) then
            Arc = 2 * F(c1, c2, pilims, r, b) - F(c1, c2, phi1, r, b) - F(c1, c2, -phi2, r, b)
            return
        else
            if (phi2 < phi1) then
                Arc = 2 * F(c1, c2, pilims, r, b) + F(c1, c2, phi2, r, b) - F(c1, c2, phi1, r, b)
                return
            else
                if (phi1 == phi2) then
                    Arc = 0.d0
                    return
                else
                    Arc = F(c1, c2, phi2, r, b) - F(c1, c2, phi1, r, b)
                    return
                end if
            end if
        end if
    end if
    
    return
    
end function

real*8 function F(c1, c2, phi, r, b)

    real*8 :: c1, c2, phi, r, b, Fc, Fq, Fl
    real*8 :: o
    real*8 :: alpha, beta, gamma, d, s, n, m, x
    real*8 :: ellipf, ellipe, ellippi, ellipf_tmp
    
    Fc = 0.5 * r * (r * phi - b * Sin(phi))
    Fq = -(r*(4*b*(2*b**2.d0 + 9*r**2.d0)*Sin(phi) &
             - 4*r*(3*(2*b**2.d0 + r**2.d0)*phi + b*r*Sin(3*phi)) &
             + r**3.d0*Sin(4*phi))) / 48.d0
    
    Fl = 0.d0
    if (b == 0.d0) then
        if (r == 1.d0) then
            Fl = phi * o3
            goto 2
        else
            Fl = phi * (1.d0 - (1.d0 - r * r) ** (1.5)) * o3
            goto 2
        end if
    else if (b == r) then
        if (r == 0.5) then
            s = phi * 0.5
            Fl = phi * o3 * 0.5 + o3 * Sin(s) &
                   * (1.d0 - Sin(s)**2.d0 * o3)
            goto 2
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
            Fl = alpha * ellipe + beta * ellipf + d
            goto 2
        end if
    else if (b + r == 1.d0) then
        s = phi * 0.5
        Fl = phi * o3 * 0.5 + Atan((2 * r - 1.d0) / Tan(s)) * o3 &
                + Atan(2 * Sqrt(r * b) * Sin(s) / (1.d0 - 2 * r)) * o3 &
                + pi *  Sqrt(1.d0 - 4 * r * b) / (2 * r - 1.d0) * o3 * 0.5 &
                + 2.d0 * o9 * Sqrt(b * r) * (2 * r * (5 * r - 2.d0) &
                - 2 * b * r * Cos(phi) - 3.d0) * Sin(s)
        goto 2
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
            s = pilims * Sign(1.d0, s) * 0.5
        else
            s = Asin(Sin(s) / Sqrt(m))
        end if

        n = 1.d0 - 1.d0 / (b - r)**2.d0
        ellipf = el1(Tan(s), Sqrt(1.d0 - m))
        o = 1.d0
        ellipe = el2(Tan(s), Sqrt(1.d0 - m), o, 1.d0 - m)
        ellippi = el3(Tan(s), Sqrt(1.d0 - m), 1.d0 - n)
        
        Fl = alpha * ellipe + beta * ellipf + gamma * ellippi + d
        goto 2
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
            s = pilims * Sign(1.d0, s) * 0.5
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

end module phot