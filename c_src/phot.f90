module phot
use iso_c_binding
use ellip

implicit none

real*8, parameter :: pi = 3.14159265358979323846, pihalf = 1.5707963267948966, twopithree = 2.0943951023931953
real*8, parameter :: o3 = 0.33333333333333333333, o9 = 0.1111111111111111111, twopi = 6.283185307179586
real*8, parameter :: pithird = 1.0471975511965976, pisixth = 0.5235987755982988

contains

subroutine compute_phis(rp, rm, bp, bm, bpm, phim1, phim2, phip1, phip2, &
    phim1_bpm, phim1_rp, phim1_rm, phim2_bpm, phim2_rp, phim2_rm, phip1_bpm, &
    phip1_rp, phip1_rm, phip2_bpm, phip2_rp, phip2_rm)

    real*8 :: rp, rm, rp2, rm2, bp, bm, bpm, bpm2, denom
    real*8, intent(out) :: phim1, phim2, phip1, phip2
    real*8, intent(out) :: phim1_bpm, phim1_rp, phim1_rm, phip1_bpm, phip1_rm, phip1_rp
    real*8, intent(out) :: phim2_bpm, phim2_rp, phim2_rm, phip2_bpm, phip2_rm, phip2_rp
    real*8 :: a, b, c, thetam, thetap, tmp, area, phip, phim
    real*8 :: thetam_bpm, thetap_bpm, thetam_rm, thetap_rm, thetam_rp, thetap_rp
    real*8 :: phim_bpm, phip_bpm
    
    bpm2 = bpm * bpm
    rp2 = rp * rp
    rm2 = rm * rm
    
    a = rm
    b = bpm
    c = rp
    if (b .gt. a) then
        tmp = b
        b = a
        a = tmp
    end if
    if (c .gt. b) then
        tmp = c
        c = b
        b = tmp
    end if
    if (b .gt. a) then
        tmp = b
        b = a
        a = tmp
    end if
    area = Sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c)))
    
    thetam = Atan2(area, (rm - rp) * (rm + rp) + bpm2)
    thetam_bpm = rm * ((bpm + rm) * (rm - bpm) - rp2) / (2 * area * bpm * rp)
    thetam_rp = rm * (bpm2 + rm2 + rp2) / (2 * area * rp2)
    thetam_rm = -rm2 / (area * rp)
    
    thetap = Atan2(area, (rp - rm) * (rp + rm) + bpm2)
    thetap_bpm = rp * ((bpm + rp) * (rp - bpm) - rm2) / (2 * area * bpm * rm)
    thetap_rp = -rp2 / (area * rm)
    thetap_rm = rp * (bpm2 + rm2 + rp2) / (2 * area * rm2)
    
    a = bm
    b = bp
    c = bpm
    if (b .gt. a) then
        tmp = b
        b = a
        a = tmp
    end if
    if (c .gt. b) then
        tmp = c
        c = b
        b = tmp
    end if
    if (b .gt. a) then
        tmp = b
        b = a
        a = tmp
    end if
    area = Sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c)))
    phim = Atan2(area, (bm - bp) * (bm + bp) + bpm2)
    phim_bpm = bm * ((bpm + bm) * (bm - bpm) - bp * bp) / (2 * area * bpm * bp)
    
    phip = Atan2(area, (bp - bm) * (bp + bm) + bpm2)
    phip_bpm = bm * ((bpm + bp) * (bp - bpm) - bm * bm) / (2 * area * bpm * bm)

    phim1 = phim + thetam
    phim1_bpm = phim_bpm + thetam_bpm
    phim1_rp = thetam_rp
    phim1_rm = thetam_rm
    
    phim2 = phim - thetam
    phim2_bpm = phim_bpm - thetam_bpm
    phim2_rp = thetam_rp
    phim2_rm = thetam_rm
    
    phip1 = phip + thetap
    phip1_bpm = phip_bpm + thetap_bpm
    phip1_rp = thetap_rp
    phip1_rm = thetap_rm
    
    phip2 = phip - thetap
    phip2_bpm = phip_bpm - thetap_bpm
    phip2_rp = thetap_rp
    phip2_rm = thetap_rm
    
    if (phim1 .gt. pi) then
        phim1 = phim1 - twopi
    end if
    if (phim2 .gt. pi) then
        phim2 = phim2 - twopi
    end if
    if (phip1 .gt. pi) then
        phip1 = phip1 - twopi
    end if
    if (phip2 .gt. pi) then
        phip2 = phip2 - twopi
    end if
end

subroutine compute_theta(rp, bp, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)

    real*8 :: a, b, c, rp, bp, theta, phi
    real*8 :: tmp, area, rp2, bp2
    real*8 :: theta_bp, theta_rp, phi_bp, phi_rp
    
    rp2 = rp * rp
    bp2 = bp * bp
        
    if (bp .gt. 1.d0) then
        a = bp
        b = 1.d0
        c = rp
        if (rp .gt. bp) then
            b = rp
            c = bp
        end if
    else
        a = 1.d0
        b = bp
        c = rp
    end if

    area = Sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c)))
    theta = Atan2(area, (rp - 1.d0) * (rp + 1.d0) + bp2)
    theta_bp = rp * ((bp + rp) * (rp - bp) - 1.d0) / (2 * area * bp)
    theta_rp = -rp2 / area
    
    phi = Atan2(area, (1.d0 - rp) * (1.d0 + rp) + bp2)
    phi_bp = ((bp + 1.d0) * (1.d0 - bp) - rp2) / (2 * area * bp * rp)
    phi_rp = (bp2 + 1.d0 + rp2) / (2 * area * rp2)
end

subroutine flux(c1, c2, rp, rm, bp2, bm2, bpm2, lc, j) bind(C, name="flux")

    integer (c_int), bind(C) :: j
    integer :: i
    real (c_double), bind(C) :: rp, rm
    real (c_double), bind(C), dimension(j) :: bp2, bm2, bpm2
    real (c_double), bind(C), intent(out), dimension(j) :: lc
    real (c_double), bind(C) :: c1, c2
    
    real*8 :: theta, thetapstar, thetamstar, thetapm, phim1, phim2, phip1, phip2, thetap, thetam
    real*8 :: phim1_bpm, phim1_rp, phim1_rm, phip1_bpm, phip1_rm, phip1_rp
    real*8 :: phim2_bpm, phim2_rp, phim2_rm, phip2_bpm, phip2_rm, phip2_rp
    real*8 :: phim, phip, cosphim, cosphip, costhetapm, costheta, cosphi1, cosphi2
    real*8 :: bp2i, bm2i, bpm2i, bpi, bmi, bpmi, rp2, rm2, obpi, obmi
    real*8 :: thetap_bp, thetap_rp, phip_bp, phip_rp, thetam_bp, thetam_rp, phim_bp, phim_rp
    real*8 :: theta_bp, theta_rp, phi_bp, phi_rp, phi
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
                        lc(i) = 1.d0 - Arc(c1, c2, -pi, pi, rm, bmi) * of0
                        goto 1
                    else
                        call compute_theta(rm, bmi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                        lc(i) = 2 * (Fstar(c1, c2, pi - phi) &
                                - F(c1, c2, theta, rm, bmi)) * of0
                        goto 1
                    end if
                end if
            else
                if (bmi .gt. rm + 1.d0) then
                    if (bpi + rp .lt. 1.d0) then
                        lc(i) = 1.d0 - 2 * Fcomplete(c1, c2, rp, bpi) * of0
                        goto 1
                    else
                        call compute_theta(rp, bpi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                        lc(i) = 2 * (Fstar(c1, c2, pi - phi) &
                              - F(c1, c2, theta, rp, bpi)) * of0
                        goto 1
                    end if
                else
                    if (bpi + rp .lt. 1.d0) then
                        if (bmi + rm .lt. 1.d0) then
                            lc(i) = 1.d0 - 2 * (Fcomplete(c1, c2, rm, bmi) &
                                  + Fcomplete(c1, c2, rp, bpi)) * of0
                            goto 1
                        else
                            call compute_theta(rm, bmi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                            lc(i) = 2 * (Fstar(c1, c2, pi - phi) &
                                  - F(c1, c2, theta, rm, bmi) &
                                  - Fcomplete(c1, c2, rp, bpi)) * of0
                            goto 1
                        end if
                    else
                        if (bmi + rm .lt. 1.d0) then
                            call compute_theta(rp, bpi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                            lc(i) = 2 * (Fstar(c1, c2, pi - phi) &
                                  - F(c1, c2, theta, rp, bpi) & 
                                  - Fcomplete(c1, c2, rm, bmi)) * of0
                            goto 1
                        else
                            call compute_theta(rp, bpi, thetap, phip, thetap_bp, thetap_rp, phip_bp, phip_rp)
                            call compute_theta(rm, bmi, thetam, phim, thetam_bp, thetam_rp, phim_bp, phim_rp)
                            phi = phip + phim
                            lc(i) = 2 * (Fstar(c1, c2, pi - phi) &
                                  - F(c1, c2, thetap, rp, bpi) & 
                                  - F(c1, c2, thetam, rm, bmi)) * of0
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
                        lc(i) = 1.d0 - 2 * Fcomplete(c1, c2, rm, bmi) * of0
                        goto 1
                    else
                        call compute_theta(rm, bmi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                        lc(i) = 2 * (Fstar(c1, c2, pi - phi) &
                              - F(c1, c2, theta, rm, bmi)) * of0
                        goto 1
                    end if
                end if
            else
                if (bmi .gt. rm + 1.d0) then
                    if (bpi + rp .lt. 1.d0) then
                        lc(i) = 1.d0 - 2 * Fcomplete(c1, c2, rp, bpi) * of0
                        goto 1
                    else
                        call compute_theta(rp,  bpi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                        lc(i) = 2 * (Fstar(c1, c2, pi - phi) &
                              - F(c1, c2, theta, rp, bpi)) * of0
                        goto 1
                    end if
                else
                    if (bpi + rp .lt. 1.d0) then
                        if (bmi + rm .lt. 1.d0) then
                            if (bpmi + rm .lt. rp) then
                                lc(i) = 1.d0 - 2 * Fcomplete(c1, c2, rp, bpi) * of0
                                goto 1
                            else
                                call compute_phis(rp, rm, bpi, bmi, bpmi, phim1, phim2, phip1, phip2, &
                                    phim1_bpm, phim1_rp, phim1_rm, phim2_bpm, phim2_rp, phim2_rm, phip1_bpm, &
                                    phip1_rp, phip1_rm, phip2_bpm, phip2_rp, phip2_rm)
                                lc(i) = 1.d0 - (Arc(c1, c2, phip1, phip2, rp, bpi) &
                                    + Arc(c1, c2, phim1, phim2, rm, bmi)) * of0
                                goto 1
                            end if
                        else 
                            if (bpmi + rm .lt. rp) then
                                goto 1
                            else
                                call compute_phis(rp, rm,  bpi, bmi, bpmi, phim1, phim2, phip1, phip2, &
                                    phim1_bpm, phim1_rp, phim1_rm, phim2_bpm, phim2_rp, phim2_rm, phip1_bpm, &
                                    phip1_rp, phip1_rm, phip2_bpm, phip2_rp, phip2_rm)
                                call compute_theta(rm, bmi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                                lc(i) = (2 * Fstar(c1, c2, pi - phi) &
                                  - Arc(c1, c2, -theta, phim2, rm, bmi) &
                                  - Arc(c1, c2, phim1, theta, rm, bmi) &
                                  - Arc(c1, c2, phip1, phip2, rp, bpi)) * of0
                                goto 1
                            end if
                        end if
                    else
                        if (bmi + rm .lt. 1.d0) then 
                            if (bpmi + rm .lt. rp) then
                                call compute_theta(rp, bpi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                                lc(i) = 2 * (Fstar(c1, c2, pi - phi) &
                                    - F(c1, c2, theta, rp, bpi)) * of0
                                goto 1
                            else
                                call compute_phis(rp, rm, bpi, bmi, bpmi, phim1, phim2, phip1, phip2, &
                                    phim1_bpm, phim1_rp, phim1_rm, phim2_bpm, phim2_rp, phim2_rm, phip1_bpm, &
                                    phip1_rp, phip1_rm, phip2_bpm, phip2_rp, phip2_rm)
                                call compute_theta(rp, bpi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                                lc(i) = (2 * Fstar(c1, c2, pi - phi) &
                                  - Arc(c1, c2, -theta, phip2, rp, bpi) &
                                  - Arc(c1, c2, phip1, theta, rp, bpi) &
                                  - Arc(c1, c2, phim1, phim2, rm, bmi)) * of0
                                goto 1
                            end if
                        else
                            obpi = 1.d0 / bpi
                            obmi = 1.d0 / bmi
                            if (bpmi + rm .lt. rp) then
                                call compute_theta(rp, bpi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                                lc(i) = 2 * (Fstar(c1, c2, pi - phi) &
                                    - F(c1, c2, theta, rp, bpi)) * of0
                                goto 1
                            else
                                costhetapm = ((bpi + bpmi) * (bpi - bpmi) + bmi * bmi) * 0.5 * obpi * obmi
                                cosphim = ((bmi - rm) * (bmi + rm) + 1.d0) * 0.5 * obmi
                                cosphip = ((bpi - rp) * (bpi + rp) + 1.d0) * 0.5 * obpi
                                thetapm = Acos(costhetapm)
                                phim = Acos(cosphim)
                                phip = Acos(cosphip)
                                if (thetapm + phim .lt. phip) then
                                        call compute_phis(rp, rm, bpi, bmi, bpmi, phim1, phim2, phip1, phip2, &
                                            phim1_bpm, phim1_rp, phim1_rm, phim2_bpm, phim2_rp, phim2_rm, phip1_bpm, &
                                            phip1_rp, phip1_rm, phip2_bpm, phip2_rp, phip2_rm)
                                        call compute_theta(rp,  bpi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                                        if (phip2 .gt. theta) then
                                            lc(i) = 2 * (Fstar(c1, c2, pi - phi) &
                                                - F(c1, c2, theta, rp, bpi)) * of0
                                            goto 1
                                        else
                                            lc(i) = (2 * Fstar(c1, c2, pi - thetap) &
                                              - Arc(c1, c2, -thetapstar, phip2, rp, bpi) &
                                              - Arc(c1, c2, phip1, thetapstar, rp, bpi) &
                                              - Arc(c1, c2, phim1, phim2, rm, bmi)) * of0
                                             goto 1
                                        end if
                                else if (thetapm + phip .lt. phim) then
                                    if ((bpi - rp) .lt. (bmi - rm)) then
                                        call compute_phis(rp, rm, bpi, bmi, bpmi, phim1, phim2, phip1, phip2, &
                                            phim1_bpm, phim1_rp, phim1_rm, phim2_bpm, phim2_rp, phim2_rm, phip1_bpm, &
                                            phip1_rp, phip1_rm, phip2_bpm, phip2_rp, phip2_rm)
                                        call compute_theta(rm, bmi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                                        lc(i) = (2 * Fstar(c1, c2, pi - phi) &
                                            - Arc(c1, c2, -theta, phim2, rm, bmi) &
                                            - Arc(c1, c2, phim1, theta, rm, bmi) & 
                                            - Arc(c1, c2, phip1, phip2, rp, bpi)) * of0
                                        goto 1
                                    else
                                        call compute_theta(rm, bmi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                                        lc(i) = 2 * (Fstar(c1, c2, pi - phi) &
                                          - F(c1, c2, theta, rm, bmi)) * of0
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
                                        call compute_theta(rp, bpi, thetap, phip, thetap_bp, thetap_rp, phip_bp, phip_rp)
                                        call compute_theta(rm, bmi, thetam, phim, thetam_bp, thetam_rp, phim_bp, phim_rp)
                                        phi = phip + phim
                                        lc(i) = 2 * (Fstar(c1, c2, pi - phi) &
                                          - F(c1, c2, thetap, rp, bpi) & 
                                          - F(c1, c2, thetam, rm, bmi)) * of0
                                        goto 1
                                    else if (d2 .lt. 1.d0) then
                                        call compute_phis(rp, rm, bpi, bmi, bpmi, phim1, phim2, phip1, phip2, &
                                            phim1_bpm, phim1_rp, phim1_rm, phim2_bpm, phim2_rp, phim2_rm, phip1_bpm, &
                                            phip1_rp, phip1_rm, phip2_bpm, phip2_rp, phip2_rm)
                                        call compute_theta(rm, bmi, thetam, phim, thetam_bp, thetam_rp, phim_bp, phim_rp)
                                        call compute_theta(rp, bpi, thetap, phip, thetap_bp, thetap_rp, phip_bp, phip_rp)
                                        phi = phim + phim
                                        lc(i) = (2 * Fstar(c1, c2, pi - phi) &
                                            - Arc(c1, c2, -thetam, -phim1, rm, bmi) & 
                                            - Arc(c1, c2, -phim2, thetam, rm, bmi) &
                                            - Arc(c1, c2, phip1, thetap, rp, bpi) &
                                            - Arc(c1, c2, -thetap, phip2, rp, bpi)) * of0
                                        goto 1
                                    else
                                        call compute_phis(rp, rm, bpi, bmi, bpmi, phim1, phim2, phip1, phip2, &
                                            phim1_bpm, phim1_rp, phim1_rm, phim2_bpm, phim2_rp, phim2_rm, phip1_bpm, &
                                            phip1_rp, phip1_rm, phip2_bpm, phip2_rp, phip2_rm)
                                        call compute_theta(rm, bmi, thetam, phim, thetam_bp, thetam_rp, phim_bp, phim_rp)
                                        call compute_theta(rp, bpi, thetap, phip, thetap_bp, thetap_rp, phip_bp, phip_rp)
                                        phi = 0.5 * (phip + phim + thetapm)
                                        lc(i) = (2 * Fstar(c1,c2, pi - phi) & 
                                            - Arc(c1, c2, -phim2, thetam, rm, bmi) & 
                                            - Arc(c1, c2, -thetap, phip2, rp, bpi)) * of0
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
    real*8 :: const, lin, quad
        
    if (phi1 < 0) then
        if (phi2 > 0) then
            Arc = F(c1, c2, phi2, r, b) + F(c1, c2, -phi1, r, b)
            return
        else
            if (phi2 < phi1) then
                Arc = 2 * Fcomplete(c1, c2, r, b) + F(c1, c2, -phi1, r, b) - F(c1, c2, -phi2, r, b)
                return
            else
                Arc = -F(c1, c2, -phi2, r, b) + F(c1, c2, -phi1, r, b)
                return
            end if
        end if
    else
        if (phi2 < 0) then
            Arc = 2 * Fcomplete(c1, c2, r, b) - F(c1, c2, phi1, r, b) - F(c1, c2, -phi2, r, b)
            return
        else
            if (phi2 < phi1) then
                Arc = 2 * Fcomplete(c1, c2, r, b) + F(c1, c2, phi2, r, b) - F(c1, c2, phi1, r, b)
                return
            else
                Arc = F(c1, c2, phi2, r, b) - F(c1, c2, phi1, r, b)
                return
            end if
        end if
    end if
    
    return
    
end function

real*8 function Fstar(c1, c2, phi)

    real*8 c1, c2, phi, Fc, Fq, Fl
    Fc = 0.5 * phi
    Fq = 0.25 * (phi + 0.5 * o3 * (Sin(2 * phi) - Sin(4 * phi)))
    Fl = phi * o3
    
    Fstar = (1.d0 - c1 - 2 * c2) * Fc + (c1 + 2 * c2) * Fl + c2 * Fq
    return
    
end function

real*8 function Fcomplete(c1, c2, r, b)

    real*8 :: c1, c2, r, b, Fc, Fq, Fl
    real*8 :: o, ome
    real*8 :: gamma, n, m, x, y, tans, sphi, br, bmr, bpr
    real*8 :: r2, b2, apb, apbo
    real*8 :: ellippi, eplusf
    
    r2 = r * r
    b2 = b * b
    bmr = b - r
    bpr = b + r
    br = b * r
    
    Fc = r2 * pihalf     
    Fq = pihalf * r2 * (b2 + 0.5 * r2) 
    
    if (-c1 .eq. 2 * c2) then
        Fl = 0.d0
        goto 3
        
    else
        if (b == r) then
            if (r == 0.5) then
            
                Fl = pisixth + 2.d0 * o9
                goto 3
                
            else
            
                m = 4 * r2
                apb = m * o9
                apbo = o9 * 2 * (1.d0 - m) * (1.d0 - 0.5 * m)
                eplusf = cel(Sqrt(1.d0 - m), o, apb, apbo)
                Fl = eplusf + pisixth
                goto 3
                
            end if
        else if (bpr == 1.d0) then
        
            y = Sqrt(br)
            Fl = pisixth + Atan(2 * y / (1.d0 - 2 * r)) * o3 &
                + pisixth *  Sqrt(1.d0 - 4 * r * b) / (2 * r - 1.d0) &
                + 2.d0 * o9 * y * (r * (10 * r - 4.d0) + 2 * br - 3.d0)
                
            goto 3
        else
        
            x = 1 / (9 * Sqrt((1 - bmr) * (1 + bmr)))
            apb = 2 * (b * (b2 - 4.d0) * r + r2 * ((6.d0 - 5 * b2) + 7 * br - 3 * r2) - 0.5 * 3.d0) * x
            apbo = (r2 * (12.d0 - 6 * r2) + b * (r * (8.d0 - 14 * r2) + b * (-2 * br - 10 * r2)) - 3.d0) * x
            gamma = 3 * x * bpr / bmr
            
            n = - 4 * r * b / (bmr * bmr)
            m = 4 * r * b / ((1 - bmr) * (1 + bmr))
            ome = Sqrt(1.d0 - m)
            
            o = 1.d0
            ellippi = cel((ome), 1.d0 - n, (o), (o)) 
            eplusf = cel((ome), (o), apb, apbo)
            Fl = eplusf + gamma * ellippi + pisixth * (1.d0 - Sign(1.d0, bmr))
            
            goto 3
        end if
    end if
    
    if (b .eq. 0.d0) then
        Fl = pithird * (1.d0 - (1.d0 - r2) ** (1.5))
        goto 3
    end if
    
3   Fcomplete = (1.d0 - c1 - 2 * c2) * Fc + (c1 + 2 * c2) * Fl + c2 * Fq
    return

end function

real*8 function F(c1, c2, phi, r, b)
    
    real*8 :: c1, c2, phi, cphi, r, b, Fc, Fq, Fl
    real*8 :: o
    real*8 :: gamma, d, s, n, m, x, y, ome, tans, sphi, br, bmr, bpr
    real*8 :: r2, b2, tanphi2, sinphi2, apb, apbo
    real*8 :: ellippi, eplusf
    
    r2 = r * r
    b2 = b * b
    bmr = b - r
    bpr = b + r
    br = b * r
    
    cphi = Cos(phi)
    sphi = Sin(phi)
    sinphi2 = Sin(phi * 0.5)
    
    if (phi .eq. pi) then
        F = Fcomplete(c1, c2, r, b)
        return
    end if
        
    Fc = 0.5 * (r2 * phi - br * sphi)
    Fq = r2 * r2 * 0.25 * (phi - sphi * cphi * (4 * o3 * cphi * cphi - 1.d0)) &
        + b * (br * 0.5 * (phi * r - b * sphi * o3) &
        + r2 * r * sphi * 0.25 * ((4 * cphi * cphi - 1.d0) * o3 - 3.d0))
        
    if (-c1 .eq. 2 * c2) then
        Fl = 0.d0
        goto 2
        
    else
        if (b == r) then
            if (r == 0.5) then
            
                Fl = phi * o3 * 0.5 + o3 * sinphi2 &
                    * (1.d0 - sinphi2 * sinphi2 * o3)
                goto 2
                
            else
            
                m = 4 * r2
                apb = m * o9
                apbo = o9 * 2 * (1.d0 - m) * (1.d0 - 0.5 * m)
                
                d = phi * o3 * 0.5 - 2 * r2 * sphi * Sqrt(1.d0 + 2 * r2 * (cphi - 1.d0)) * o9
                eplusf = el2(Tan(phi * 0.5), Sqrt(1.d0 - m), apb, apbo)
                Fl = eplusf + d
                goto 2
                
            end if
        else if (bpr == 1.d0) then
        
            y = Sqrt(br)
            Fl = phi * o3 * 0.5 + Atan((2 * r - 1.d0) / Tan(phi * 0.5)) * o3 &
                + Atan(2 * y * sinphi2 / (1.d0 - 2 * r)) * o3 &
                + pisixth *  Sqrt(1.d0 - 4 * r * b) / (2 * r - 1.d0) &
                + 2.d0 * o9 * y * (2 * r * (5 * r - 2.d0) &
                - 2 * br * cphi - 3.d0) * sinphi2
            goto 2
            
        else if (bpr .gt. 1.d0) then
        
            y = 1.d0 / (18 * Sqrt(br))
            apb = (r2 * (12.d0 - 6 * r2) + b * (b * (2 * b * r - 10 * r2) + r * (14 * r2 - 8.d0)) - 3.d0) * y
            apbo = (1.d0 + b2 * (b2 - 2 * r2 - 5.d0) + r2 * (1.d0 + r2)) * y
            gamma = bpr * y * 3 / bmr
            
            m = (1 - bmr) * (1 + bmr) / (4 * r * b)
            ome = Sqrt((b + r - 1.d0) * (b + r + 1.d0) / (4 * b * r))
            
            if (abs(sinphi2 / Sqrt(m) - 1.D0) .lt. 1.D-12) then
            
                d = phi * o3 * 0.5 - Atan2(bpr * sinphi2, bmr * Cos(phi * 0.5)) * o3
                n = 1.d0 - 1.d0 / (bmr * bmr)
                o = 1.d0
                ellippi = cel((ome), 1.d0 - n, (o), (o))
                eplusf = cel((ome), (o), apb, apbo)
                
            else
            
                d = o3 * (phi * 0.5 - Atan2(bpr * sinphi2, bmr * Cos(phi * 0.5))) &
                    - (2 * br * o9) * sphi * Sqrt(1.d0 - (b2 + r2) + 2 * br * cphi)
                tans = 1.d0 / Sqrt(m / (sinphi2 * sinphi2) - 1.d0)
                n = 1.d0 - 1.d0 / (bmr * bmr)
                ellippi = el3((tans), (ome), 1.d0 - n)
                eplusf = el2((tans), (ome), apb, apbo)
                
            end if
            
            Fl = eplusf + gamma * ellippi + d
            goto 2
            
        else
        
            x = 1 / (9 * Sqrt((1 - bmr) * (1 + bmr)))
            apb = 2 * (b * (b2 - 4.d0) * r + r2 * ((6.d0 - 5 * b2) + 7 * br - 3 * r2) - 0.5 * 3.d0) *x
            apbo = (r2 * (12.d0 - 6 * r2) + b * (r * (8.d0 - 14 * r2) + b * (-2 * br - 10 * r2)) - 3.d0) *x
            gamma = 3 * x * bpr / bmr
            
            n = - 4 * r * b / (bmr * bmr)
            m = 4 * r * b / ((1 - bmr) * (1 + bmr))
            ome = Sqrt(1.d0 - m)
            
            tanphi2 = Tan(phi * 0.5)
            d = o3 * (phi * 0.5 - Atan((bpr / bmr) * tanphi2)) &
                - 2 * br * o9 * sphi * Sqrt(1.d0 - (b2 + r2) + 2 * br * cphi)
            ellippi = el3((tanphi2), (ome), 1.d0 - n)
            eplusf = el2((tanphi2), (ome), apb, apbo)
            
            Fl = eplusf + gamma * ellippi + d
            goto 2
            
        end if
        
        if (b .eq. 0.d0) then
            Fl = phi * (1.d0 - (1.d0 - r2) ** (1.5)) * o3
            goto 2
        end if
        
    end if
    
2   F = (1.d0 - c1 - 2 * c2) * Fc + (c1 + 2 * c2) * Fl + c2 * Fq

    return
end function

end module phot