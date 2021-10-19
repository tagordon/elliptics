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
    real*8 :: phim_bpm, phip_bpm, x, y
    
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
    thetam_bpm = ((bpm + rm) * (rm - bpm) - rp2) / (area * bpm)
    thetam_rp = 2 * rp / area
    thetam_rm = ((bpm - rm) * (bpm + rm) - rp2) / (area * rm)
    
    thetap = Atan2(area, (rp - rm) * (rp + rm) + bpm2)
    thetap_bpm = ((rp - rm) * (rp + rm) - bpm2) / (bpm * area)
    thetap_rp = ((bpm - rp) * (bpm + rp) - rm2) / (rp * area)
    thetap_rm = 2 * rm / area
    
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
    phim_bpm = ((bm - bp) * (bm + bp) - bpm2) / (bpm * area)
    
    phip = Atan2(area, (bp - bm) * (bp + bm) + bpm2)
    phip_bpm = ((bp - bm) * (bp + bm) - bpm2) / (bpm * area)

    phim1 = phim + thetam
    phim1_bpm = phim_bpm + thetam_bpm
    phim1_rp = thetam_rp
    phim1_rm = thetam_rm
    
    phim2 = phim - thetam
    phim2_bpm = phim_bpm - thetam_bpm
    phim2_rp = -thetam_rp
    phim2_rm = -thetam_rm
    
    phip1 = phip + thetap
    phip1_bpm = phip_bpm + thetap_bpm
    phip1_rp = thetap_rp
    phip1_rm = thetap_rm
    
    phip2 = phip - thetap
    phip2_bpm = phip_bpm - thetap_bpm
    phip2_rp = -thetap_rp
    phip2_rm = -thetap_rm
    
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
    real*8 :: tmp, area, rp2, bp2, x
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
    theta_bp = ((rp + bp) * (rp - bp) - 1.d0) / (bp * area)
    theta_rp = ((bp + rp) * (bp - rp) - 1.d0) / (rp * area)
    
    phi = Atan2(area, (1.d0 - rp) * (1.d0 + rp) + bp2)
    phi_bp = -((bp + 1.d0) * (bp - 1.d0) + rp2) / (bp * area)
    phi_rp = 2 * rp / area
end

subroutine flux(c1, c2, rp, rm, bp2, bm2, bpm2, lc, j) bind(C, name="flux")

    integer (c_int), bind(C) :: j
    integer :: i
    real (c_double), bind(C) :: rp, rm
    real (c_double), bind(C), dimension(j) :: bp2, bm2, bpm2
    real (c_double), bind(C), intent(out), dimension(8, j) :: lc
    real*8, dimension(8) :: f0
    real (c_double), bind(C) :: c1, c2
    
    real*8 :: theta, thetapm, phim1, phim2, phip1, phip2, thetap, thetam
    real*8 :: phim1_bpm, phim1_rp, phim1_rm, phip1_bpm, phip1_rm, phip1_rp
    real*8 :: phim2_bpm, phim2_rp, phim2_rm, phip2_bpm, phip2_rm, phip2_rp
    real*8 :: phim, phip, cosphim, cosphip, costhetapm, costheta, cosphi1, cosphi2
    real*8 :: bp2i, bm2i, bpm2i, bpi, bmi, bpmi, rp2, rm2, obpi, obmi
    real*8 :: thetap_bp, thetap_rp, phip_bp, phip_rp, thetam_bp, thetam_rp, phim_bp, phim_rp
    real*8 :: theta_bp, theta_rp, phi_bp, phi_rp, phi, thetam_bm, thetam_rm, phi_bm, phi_rm
    real*8 :: d1, d2, area, phim_bm, phim_rm, theta_bm, theta_rm, theta_bpm, thetapm_bpm
    real*8 :: of0, a, b, c, tmp
    
    real*8 :: bp(j), bm(j), bpm(j)
    bp = Sqrt(bp2)
    bm = Sqrt(bm2)
    bpm = Sqrt(bpm2)
    
    rp2 = rp * rp
    rm2 = rm * rm
    
    f0(1) = (1.d0 - c1 - 2 * c2) * pi + (c1 + 2 * c2) * twopithree + c2 * pihalf
    f0(2) = 0.d0
    f0(3) = 0.d0
    f0(4) = 0.d0
    f0(5) = 0.d0
    f0(6) = 0.d0
    f0(7) = -pi + twopithree 
    f0(8) = -2 * pi + 2 * twopithree + pihalf
    
    of0 = 1.d0 / f0(1)
    lc = 0.d0
    
    do i=1,j,1
    
        bp2i = bp2(i)
        bm2i = bm2(i)
        bpm2i = bpm2(i)
        bpi = bp(i)
        bmi = bm(i)
        bpmi = bpm(i)
        
        if ((bpi .gt. rp + 1.d0) .AND. (bmi .gt. rm + 1.d0)) then
            lc(:, i) = f0 * of0
            goto 1
        else if (bpmi .gt. rp + rm) then
            if (bpi .gt. rp + 1.d0) then
                if (bmi .gt. rm + 1.d0) then
                    lc(:, i) = f0 * of0
                    goto 1
                else
                    if (bmi + rm .lt. 1.d0) then
                        lc(:, i) = (f0 - 2 * Fcomplete(c1, c2, rm, bmi, .FALSE.)) * of0
                        goto 1
                    else
                        call compute_theta(rm, bmi, theta, phi, theta_bm, theta_rm, phi_bm, phi_rm)
                        lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, 0.d0, 0.d0, -phi_bm, -phi_rm, 0.d0) &
                                - F(c1, c2, theta, rm, bmi, 0.d0, 0.d0, theta_bm, theta_rm, 0.d0, .FALSE.)) * of0
                        goto 1
                    end if
                end if
            else
                if (bmi .gt. rm + 1.d0) then
                    if (bpi + rp .lt. 1.d0) then
                        lc(:, i) = (f0 - 2 * Fcomplete(c1, c2, rp, bpi, .TRUE.)) * of0
                        goto 1
                    else
                        call compute_theta(rp, bpi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                        lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, -phi_bp, -phi_rp, 0.d0, 0.d0, 0.d0) &
                              - F(c1, c2, theta, rp, bpi, theta_bp, theta_rp, 0.d0, 0.d0, 0.d0, .TRUE.)) * of0
                        goto 1
                    end if
                else
                    if (bpi + rp .lt. 1.d0) then
                        if (bmi + rm .lt. 1.d0) then
                            lc(:, i) = (f0 - 2 * (Fcomplete(c1, c2, rm, bmi, .FALSE.) &
                                  + Fcomplete(c1, c2, rp, bpi, .TRUE.))) * of0
                            goto 1
                        else
                            call compute_theta(rm, bmi, theta, phi, theta_bm, theta_rm, phi_bm, phi_rm)
                            lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, 0.d0, 0.d0, -phi_bm, -phi_rm, 0.d0) &
                                  - F(c1, c2, theta, rm, bmi, 0.d0, 0.d0, theta_bm, theta_rm, 0.d0, .FALSE.) &
                                  - Fcomplete(c1, c2, rp, bpi, .TRUE.)) * of0
                            goto 1
                        end if
                    else
                        if (bmi + rm .lt. 1.d0) then
                            call compute_theta(rp, bpi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                            lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, -phi_bp, -phi_rp, 0.d0, 0.d0, 0.d0) &
                                  - F(c1, c2, theta, rp, bpi, theta_bp, theta_rp, 0.d0, 0.d0, 0.d0, .TRUE.) & 
                                  - Fcomplete(c1, c2, rm, bmi, .FALSE.)) * of0
                            goto 1
                        else
                            call compute_theta(rp, bpi, thetap, phip, thetap_bp, thetap_rp, phip_bp, phip_rp)
                            call compute_theta(rm, bmi, thetam, phim, thetam_bm, thetam_rm, phim_bm, phim_rm)
                            phi = phip + phim
                            !phi_bp = phip_bp + phim_bp
                            !phi_rp = phip_rp + phim_rp
                            lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, -phip_bp, -phip_rp, -phim_bm, -phim_rm, 0.d0) &
                                  - F(c1, c2, thetap, rp, bpi, thetap_bp, thetap_rp, 0.d0, 0.d0, 0.d0, .TRUE.) & 
                                  - F(c1, c2, thetam, rm, bmi, 0.d0, 0.d0, thetam_bm, thetam_rm, 0.d0, .FALSE.)) * of0
                            !lc(:, i) = 1.1d0
                            goto 1
                        end if
                    end if
                end if
            end if
        else
            if (bpi .gt. rp + 1.d0) then
                if (bmi .gt. rm + 1.d0) then
                    lc(:, i) = f0
                    goto 1
                else
                    call compute_theta(rm, bmi, theta, phi, theta_bm, theta_rm, phi_bm, phi_rm)
                    lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, 0.d0, 0.d0, -phi_bm, -phi_rm, 0.d0) &
                            - F(c1, c2, theta, rm, bmi, 0.d0, 0.d0, theta_bm, theta_rm, 0.d0, .FALSE.)) * of0
                    goto 1
                end if
            else
                if (bmi .gt. rm + 1.d0) then
                    if (bpi + rp .lt. 1.d0) then
                        lc(:, i) = (f0 - 2 * Fcomplete(c1, c2, rp, bpi, .TRUE.)) * of0
                        goto 1
                    else
                        call compute_theta(rp,  bpi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                        lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, -phi_bp, -phi_rp, 0.d0, 0.d0, 0.d0) &
                              - F(c1, c2, theta, rp, bpi, theta_bp, theta_rp, 0.d0, 0.d0, 0.d0, .TRUE.)) * of0
                        goto 1
                    end if
                else
                    if (bpi + rp .lt. 1.d0) then
                        if (bmi + rm .lt. 1.d0) then
                            if (bpmi + rm .lt. rp) then
                                lc(:, i) = (f0 - 2 * Fcomplete(c1, c2, rp, bpi, .TRUE.)) * of0
                                goto 1
                            else
                                call compute_phis(rp, rm, bpi, bmi, bpmi, phim1, phim2, phip1, phip2, &
                                    phim1_bpm, phim1_rp, phim1_rm, phim2_bpm, phim2_rp, phim2_rm, phip1_bpm, &
                                    phip1_rp, phip1_rm, phip2_bpm, phip2_rp, phip2_rm)
                                lc(:, i) = (f0 - (Arc(c1, c2, phip1, phip2, rp, bpi, &
                                                    0.d0, phip1_rp, 0.d0, phip1_rm, phip1_bpm, &
                                                    0.d0, phip2_rp, 0.d0, phip2_rm, phip2_bpm, .TRUE.) &
                                                + Arc(c1, c2, phim1, phim2, rm, bmi, &
                                                      0.d0, phim1_rp, 0.d0, phim1_rm, phim1_bpm, &
                                                      0.d0, phim2_rp, 0.d0, phim2_rm, phim2_bpm, .False.))) * of0
                                goto 1
                            end if
                        else 
                            if (bpmi + rm .lt. rp) then
                                goto 1
                            else
                                call compute_phis(rp, rm,  bpi, bmi, bpmi, phim1, phim2, phip1, phip2, &
                                    phim1_bpm, phim1_rp, phim1_rm, phim2_bpm, phim2_rp, phim2_rm, phip1_bpm, &
                                    phip1_rp, phip1_rm, phip2_bpm, phip2_rp, phip2_rm)
                                call compute_theta(rm, bmi, theta, phi, theta_bm, theta_rm, phi_bm, phi_rm)
                                lc(:, i) = (2 * Fstar(c1, c2, pi - phi, 0.d0, 0.d0, -phi_bm, -phi_rm, 0.d0) &
                                  - Arc(c1, c2, -theta, phim2, rm, bmi, &
                                        0.d0, 0.d0, -theta_bm, -theta_rm, 0.d0, &
                                        0.d0, phim2_rp, 0.d0, phim2_rm, phim2_bpm, .FALSE.) &
                                  - Arc(c1, c2, phim1, theta, rm, bmi, &
                                        0.d0, phim1_rp, 0.d0, phim1_rm, phim1_bpm, &
                                        0.d0, 0.d0, theta_bm, theta_rm, 0.d0, .FALSE.) &
                                  - Arc(c1, c2, phip1, phip2, rp, bpi, &
                                        0.d0, phip1_rp, 0.d0, phip1_rm, phip1_bpm, &
                                        0.d0, phip2_rp, 0.d0, phip2_rm, phip2_bpm, .TRUE.)) * of0
                                goto 1
                            end if
                        end if
                    else
                        if (bmi + rm .lt. 1.d0) then 
                            if (bpmi + rm .lt. rp) then
                                call compute_theta(rp, bpi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                                lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, -phi_bp, -phi_rp, 0.d0, 0.d0, 0.d0) &
                                    - F(c1, c2, theta, rp, bpi, theta_bp, theta_rp, 0.d0, 0.d0, 0.d0, .TRUE.)) * of0
                                goto 1
                            else
                                call compute_phis(rp, rm, bpi, bmi, bpmi, phim1, phim2, phip1, phip2, &
                                    phim1_bpm, phim1_rp, phim1_rm, phim2_bpm, phim2_rp, phim2_rm, phip1_bpm, &
                                    phip1_rp, phip1_rm, phip2_bpm, phip2_rp, phip2_rm)
                                call compute_theta(rp, bpi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                                lc(:, i) = (2 * Fstar(c1, c2, pi - phi, -phi_bp, -phi_rp, 0.d0, 0.d0, 0.d0) &
                                  - Arc(c1, c2, -theta, phip2, rp, bpi, &
                                        -theta_bp, -theta_rp, 0.d0, 0.d0, 0.d0, &
                                        0.d0, phip2_rp, 0.d0, phip2_rm, phip2_bpm, .TRUE.) &
                                  - Arc(c1, c2, phip1, theta, rp, bpi, &
                                        0.d0, phip1_rp, 0.d0, phip1_rm, phip1_bpm, &
                                        theta_bp, theta_rp, 0.d0, 0.d0, 0.d0, .TRUE.) &
                                  - Arc(c1, c2, phim1, phim2, rm, bmi, & 
                                        0.d0, phim1_rp, 0.d0, phim1_rm, phim1_bpm, &
                                        0.d0, phim2_rp, 0.d0, phim2_rm, phim2_bpm, .FALSE.)) * of0
                                goto 1
                            end if
                        else
                            obpi = 1.d0 / bpi
                            obmi = 1.d0 / bmi
                            if (bpmi + rm .lt. rp) then
                                call compute_theta(rp, bpi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                                lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, -phi_bp, -phi_rp, 0.d0, 0.d0, 0.d0) &
                                    - F(c1, c2, theta, rp, bpi, theta_bp, theta_rp, 0.d0, 0.d0, 0.d0, .TRUE.)) * of0
                                goto 1
                            else
                                call compute_theta(rp,  bpi, theta, phip, theta_bp, theta_rp, phip_bp, phip_rp)
                                call compute_theta(rm,  bmi, theta, phim, theta_bm, theta_rm, phim_bm, phim_rm)
                                
                                a = bmi
                                b = bpi
                                c = bpmi
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
                                thetapm = Atan2(area, (bmi - bpmi) * (bmi + bpmi) + bpi * bpi)
                                thetapm_bpm = 2 * bpmi / area
                                
                                if (thetapm + phim .lt. phip) then
                                        call compute_phis(rp, rm, bpi, bmi, bpmi, phim1, phim2, phip1, phip2, &
                                            phim1_bpm, phim1_rp, phim1_rm, phim2_bpm, phim2_rp, phim2_rm, phip1_bpm, &
                                            phip1_rp, phip1_rm, phip2_bpm, phip2_rp, phip2_rm)
                                        call compute_theta(rp,  bpi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                                        if (phip2 .gt. theta) then
                                            lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, -phi_bp, -phi_rp, 0.d0, 0.d0, 0.d0) &
                                                - F(c1, c2, theta, rp, bpi, theta_bp, theta_rp, 0.d0, 0.d0, 0.d0, .TRUE.)) * of0
                                            goto 1
                                        else
                                            lc(:, i) = (2 * Fstar(c1, c2, pi - phi, -phi_bp, -phi_rp, 0.d0, 0.d0, 0.d0) &
                                              - Arc(c1, c2, -theta, phip2, rp, bpi, &
                                                    -theta_bp, -theta_rp, 0.d0, 0.d0, 0.d0, &
                                                    0.d0, phip2_rp, 0.d0, phip2_rm, phip2_bpm, .TRUE.) &
                                              - Arc(c1, c2, phip1, theta, rp, bpi, &
                                                    0.d0, phip1_rp, 0.d0, phip1_rm, phip1_bpm, &
                                                    theta_bp, theta_rp, 0.d0, 0.d0, 0.d0, .TRUE.) &
                                              - Arc(c1, c2, phim1, phim2, rm, bmi, &
                                                    0.d0, phim1_rp, 0.d0, phim1_rm, phim1_bpm, &
                                                    0.d0, phim2_rp, 0.d0, phim2_rm, phim2_bpm, .FALSE.)) * of0
                                             goto 1
                                        end if
                                else if (thetapm + phip .lt. phim) then
                                    if ((bpi - rp) .lt. (bmi - rm)) then
                                        call compute_phis(rp, rm, bpi, bmi, bpmi, phim1, phim2, phip1, phip2, &
                                            phim1_bpm, phim1_rp, phim1_rm, phim2_bpm, phim2_rp, phim2_rm, phip1_bpm, &
                                            phip1_rp, phip1_rm, phip2_bpm, phip2_rp, phip2_rm)
                                        call compute_theta(rm, bmi, theta, phi, theta_bm, theta_rm, phi_bm, phi_rm)
                                        lc(:, i) = (2 * Fstar(c1, c2, pi - phi, 0.d0, 0.d0, -phi_bm, -phi_rm, 0.d0) &
                                            - Arc(c1, c2, -theta, phim2, rm, bmi, &
                                                  0.d0, 0.d0, -theta_bm, -theta_rm, 0.d0, &
                                                  0.d0, phim2_rp, 0.d0, phim2_rm, phim2_bpm, .FALSE.) &
                                            - Arc(c1, c2, phim1, theta, rm, bmi, &
                                                  0.d0, phim1_rp, 0.d0, phim1_rm, phim1_bpm, &
                                                  0.d0, 0.d0, theta_bm, theta_rm, 0.d0, .FALSE.) & 
                                            - Arc(c1, c2, phip1, phip2, rp, bpi, &
                                                  0.d0, phip1_rp, 0.d0, phip1_rm, phip1_bpm, &
                                                  0.d0, phip2_rp, 0.d0, phip2_rm, phip2_bpm, .TRUE.)) * of0
                                        goto 1
                                    else
                                        call compute_theta(rm, bmi, theta, phi, theta_bm, theta_rm, phi_bm, phi_rm)
                                        lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, 0.d0, 0.d0, -phi_bm, -phi_rm, 0.d0) &
                                          - F(c1, c2, theta, rm, bmi, 0.d0, 0.d0, theta_bm, theta_rm, 0.d0, .FALSE.)) * of0
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
                                        call compute_theta(rm, bmi, thetam, phim, thetam_bm, thetam_rm, phim_bm, phim_rm)
                                        phi = phip + phim
                                        lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, -phip_bp, -phip_rp, -phim_bm, -phim_rm, 0.d0) &
                                          - F(c1, c2, thetap, rp, bpi, theta_bp, thetap_rp, 0.d0, 0.d0, 0.d0, .TRUE.) & 
                                          - F(c1, c2, thetam, rm, bmi, 0.d0, 0.d0, thetam_rm, thetam_rp, 0.d0, .FALSE.)) * of0
                                        goto 1
                                    else if (d2 .lt. 1.d0) then
                                        call compute_phis(rp, rm, bpi, bmi, bpmi, phim1, phim2, phip1, phip2, &
                                            phim1_bpm, phim1_rp, phim1_rm, phim2_bpm, phim2_rp, phim2_rm, phip1_bpm, &
                                            phip1_rp, phip1_rm, phip2_bpm, phip2_rp, phip2_rm)
                                        call compute_theta(rm, bmi, thetam, phim, thetam_bm, thetam_rm, phim_bm, phim_rm)
                                        call compute_theta(rp, bpi, thetap, phip, thetap_bp, thetap_rp, phip_bp, phip_rp)
                                        phi = phip + phim
                                        lc(:, i) = (2 * Fstar(c1, c2, pi - phi, -phip_bp, -phip_rp, -phim_bm, -phim_rm, 0.d0) &
                                            - Arc(c1, c2, -thetam, -phim1, rm, bmi, &
                                                  0.d0, 0.d0, -thetam_bm, -thetam_rm, 0.d0, &
                                                  0.d0, -phim1_rp, 0.d0, -phim1_rm, -phim1_bpm, .FALSE.) & 
                                            - Arc(c1, c2, -phim2, thetam, rm, bmi, &
                                                  0.d0, -phim2_rp, 0.d0, -phim2_rm, -phim2_bpm, &
                                                  0.d0, 0.d0, -theta_bm, -theta_rm, 0.d0, .FALSE.) &
                                            - Arc(c1, c2, phip1, thetap, rp, bpi, &
                                                  0.d0, phip1_rp, 0.d0, phip1_rm, phip1_bpm, &
                                                  thetap_bp, thetap_rp, 0.d0, 0.d0, 0.d0, .TRUE.) &
                                            - Arc(c1, c2, -thetap, phip2, rp, bpi, &
                                                  -thetap_bp, -thetap_rp, 0.d0, 0.d0, 0.d0, &
                                                  0.d0, phip2_rp, 0.d0, phip2_rm, phip2_bpm, .TRUE.)) * of0
                                        goto 1
                                    else
                                        ! this is a place to check when things are broken later
                                        call compute_phis(rp, rm, bpi, bmi, bpmi, phim1, phim2, phip1, phip2, &
                                            phim1_bpm, phim1_rp, phim1_rm, phim2_bpm, phim2_rp, phim2_rm, phip1_bpm, &
                                            phip1_rp, phip1_rm, phip2_bpm, phip2_rp, phip2_rm)
                                        call compute_theta(rm, bmi, thetam, phim, thetam_bm, thetam_rm, phim_bm, phim_rm)
                                        call compute_theta(rp, bpi, thetap, phip, thetap_bp, thetap_rp, phip_bp, phip_rp)
                                        phi = 0.5 * (phip + phim + thetapm)
                                        lc(:, i) = (2 * Fstar(c1, c2, pi - phi, -0.5 * phip_bp, -0.5 * phip_rp, &
                                                           -0.5 * phim_bm, -0.5 * phim_rm, -0.5 * thetapm_bpm) & 
                                            - Arc(c1, c2, -phim2, thetam, rm, bmi, &
                                                  0.d0, -phim2_rp, 0.d0, -phim2_rm, -phim2_bpm, &
                                                  0.d0, 0.d0, thetam_bm, thetam_rm, 0.d0, .FALSE.) & 
                                            - Arc(c1, c2, -thetap, phip2, rp, bpi, &
                                                  -thetap_bp, -thetap_rp, 0.d0, 0.d0, 0.d0, &
                                                  0.d0, phip2_rp, 0.d0, phip2_rm, phip2_bpm, .TRUE.)) * of0
                                        goto 1
                                    end if
                                end if
                            end if
                        end if
                    end if
                end if
            end if  
1        end if
        lc(:, i) = lc(:, i) - f0 * of0
    end do
    return
    
end

function Arc(c1, c2, phi1, phi2, r, b, phi1_bp, phi1_rp, phi1_bm, phi1_rm, phi1_bpm, &
            phi2_bp, phi2_rp, phi2_bm, phi2_rm, phi2_bpm, pflag)
                    
    real*8, dimension(8) :: Arc

    logical :: pflag
    real*8 :: phi1, phi2, r, b, c1, c2
    real*8 :: phi1_bp, phi1_rp, phi1_bm, phi1_rm, phi1_bpm
    real*8 :: phi2_bp, phi2_rp, phi2_bm, phi2_rm, phi2_bpm
    real*8 :: const, lin, quad
        
    if (phi1 < 0) then
        if (phi2 > 0) then
            Arc = F(c1, c2, phi2, r, b, phi2_bp, phi2_rp, phi2_bm, phi2_rm, phi2_bpm, pflag) &
                + F(c1, c2, -phi1, r, b, -phi1_bp, -phi1_rp, -phi1_bm, -phi1_rm, -phi1_bpm, pflag)
            return
        else
            if (phi2 < phi1) then
                Arc = 2 * Fcomplete(c1, c2, r, b, pflag) &
                    + F(c1, c2, -phi1, r, b, -phi1_bp, -phi1_rp, -phi1_bm, -phi1_rm, -phi1_bpm, pflag) &
                    - F(c1, c2, -phi2, r, b, -phi2_bp, -phi2_rp, -phi2_bm, -phi2_rm, -phi2_bpm, pflag)
                return
            else
                Arc = - F(c1, c2, -phi2, r, b, -phi2_bp, -phi2_rp, -phi2_bm, -phi2_rm, -phi2_bpm, pflag) &
                    + F(c1, c2, -phi1, r, b, -phi1_bp, -phi1_rp, -phi1_bm, -phi1_rm, -phi1_bpm, pflag)
                return
            end if
        end if
    else
        if (phi2 < 0) then
            Arc = 2 * Fcomplete(c1, c2, r, b, pflag) &
                - F(c1, c2, phi1, r, b, phi1_bp, phi1_rp, phi1_bm, phi1_rm, phi1_bpm, pflag) &
                - F(c1, c2, -phi2, r, b, -phi2_bp, -phi2_rp, -phi2_bm, -phi2_rm, -phi2_bpm, pflag)
            return
        else
            if (phi2 < phi1) then
                Arc = 2 * Fcomplete(c1, c2, r, b, pflag) &
                    + F(c1, c2, phi2, r, b, phi2_bp, phi2_rp, phi2_bm, phi2_rm, phi2_bpm, pflag) &
                    - F(c1, c2, phi1, r, b, phi1_bp, phi1_rp, phi1_bm, phi1_rm, phi1_bpm, pflag)
                return
            else
                Arc = F(c1, c2, phi2, r, b, phi2_bp, phi2_rp, phi2_bm, phi2_rm, phi2_bpm, pflag) &
                    - F(c1, c2, phi1, r, b, phi1_bp, phi1_rp, phi1_bm, phi1_rm, phi1_bpm, pflag)
                return
            end if
        end if
    end if
    
    return
    
end function

function Fstar(c1, c2, phi, phi_bp, phi_rp, phi_bm, phi_rm, phi_bpm)

    real*8, dimension(8) :: Fstar

    real*8 :: c1, c2, phi, Fc, Fq, Fl
    real*8 :: Fc_bp, Fc_rp, Fc_bm, Fc_rm, Fc_bpm
    real*8 :: Fq_bp, Fq_rp, Fq_bm, Fq_rm, Fq_bpm
    real*8 :: Fl_bp, Fl_rp, Fl_bm, Fl_rm, Fl_bpm
    real*8 :: Fc_phi, Fq_phi, Fl_phi
    real*8 :: phi_bp, phi_rp, phi_bm, phi_rm, phi_bpm
    
    Fc = 0.5 * phi
    Fc_phi = 0.5
    Fc_bp = Fc_phi * phi_bp
    Fc_rp = Fc_phi * phi_rp
    Fc_bm = Fc_phi * phi_bm
    Fc_rm = Fc_phi * phi_rm
    Fc_bpm = Fc_phi * phi_bpm
    
    Fq = 0.25 * (phi + 0.5 * o3 * (Sin(2 * phi) - Sin(4 * phi)))
    Fq_phi = 0.5 * o3 * Cos(phi)**2 * (5.d0 - 4 * Cos(2 * phi))
    Fq_bp = Fq_phi * phi_bp
    Fq_rp = Fq_phi * phi_rp
    Fq_bm = Fq_phi * phi_bm
    Fq_rm = Fq_phi * phi_rm
    Fq_bpm = Fq_phi * phi_bpm
    
    Fl = phi * o3
    Fl_phi = o3
    Fl_bp = Fl_phi * phi_bp
    Fl_rp = Fl_phi * phi_rp
    Fl_bm = Fl_phi * phi_bm
    Fl_rm = Fl_phi * phi_rm
    Fl_bpm = Fl_phi * phi_bpm
    
    Fstar(1) = (1.d0 - c1 - 2 * c2) * Fc + (c1 + 2 * c2) * Fl + c2 * Fq
    Fstar(2) = (1.d0 - c1 - 2 * c2) * Fc_bp + (c1 + 2 * c2) * Fl_bp + c2 * Fq_bp
    Fstar(3) = (1.d0 - c1 - 2 * c2) * Fc_rp + (c1 + 2 * c2) * Fl_rp + c2 * Fq_rp
    Fstar(4) = (1.d0 - c1 - 2 * c2) * Fc_bm + (c1 + 2 * c2) * Fl_bm + c2 * Fq_bm
    Fstar(5) = (1.d0 - c1 - 2 * c2) * Fc_rm + (c1 + 2 * c2) * Fl_rm + c2 * Fq_rm
    Fstar(6) = (1.d0 - c1 - 2 * c2) * Fc_bpm + (c1 + 2 * c2) * Fl_bpm + c2 * Fq_bpm
    Fstar(7) = - Fc + Fl 
    Fstar(8) = -2 * Fc + 2 * Fl + Fq
    return
    
end function

function Fcomplete(c1, c2, r, b, pflag)

    logical :: pflag
    real*8, dimension(8) :: Fcomplete
    
    real*8 :: c1, c2, r, b, Fc, Fq, Fl
    real*8 :: Fc_r, Fq_r, Fl_r, Fc_b, Fq_b, Fl_b
    real*8 :: o, ome
    real*8 :: gamma, n, m, x, y, tans, sphi, br, bmr, bpr
    real*8 :: r2, b2, apb, apbo
    real*8 :: ellippi, eplusf
    real*8 :: ellippi_r, eplusf_r, alpha, alpha_r, beta, beta_r
    real*8 :: u, v, p, q, m_r, n_r, gamma_r
    
    r2 = r * r
    b2 = b * b
    bmr = b - r
    bpr = b + r
    br = b * r
    
    Fc = r2 * pihalf  
    Fc_r = r * pi
    Fc_b = 0.d0
    
    Fq = pihalf * r2 * (b2 + 0.5 * r2) 
    Fq_r = pi * r * (b2 + r2)
    Fq_b = pi * r2 * b
    
    if (-c1 .eq. 2 * c2) then
        Fl = 0.d0
        Fl_r = 0.d0
        Fl_b = 0.d0
        goto 3
        
    else
        if (b == r) then
            if (r == 0.5) then
            
                Fl = pisixth + 2.d0 * o9
                Fl_r = 0.d0
                Fl_b = 0.d0
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
                + pisixth *  Sqrt(1.d0 - 4 * br) / (2 * r - 1.d0) &
                + 2.d0 * o9 * y * (r * (10 * r - 4.d0) + 2 * br - 3.d0)
                
            Fl_r = 0.5 * o9 * (8 * y * (10 * r + b - 2.d0) &
                + 6 * pi * b / ((1.d0 - 2 * r) * Sqrt(1.d0 - 4 * br)) &
                - 6 * pi * Sqrt(1.d0 - 4 * br) / (1 - 2 * r)**2.d0 &
                + pi * (b + 2 * br) / (y * (1.d0 - 4 * r * (b + r - 1))) &
                + 2 * b * (2 * r * (5 * r + b - 2) - 3.d0) / y)
                
            Fl_b = 0.5 * o9 * r * (8 * y + 6 * pi / ((1.d0 - 2 * r) * Sqrt(1.d0 - 4 * r)) &
                + (pi - 2 * pi * r) / (y * (1.d0 + 4 * r * (b + r - 1))) &
                + (4 * r * (5 * r + b - 2.d0) - 6.d0) / y)
                
            goto 3
        else
        
            x = Sqrt((1.d0 - bmr) * (1.d0 + bmr))
            
            alpha = (7 * r2 + b2 - 4.d0) * x * o9
            alpha_r = (b * (b2 - 15 * br) + 3 * r * (6.d0 - 7 * r2) + b * (-4.d0 + 35 * r2)) &
                    / (9 * x)
            
            beta = (r2 * r2 + b2 * b2 + r2 - b2 * (5.d0 + 2 * r2) + 1.d0) / (9.d0 * x)
            beta_r = (2 * (1.d0 - bmr) * (1.d0 + bmr) * r * (1.d0 - 2 * b2 + 2 * r2) &
                    - bmr * (1.d0 + b2 * b2 + r2 + r2 * r2 &
                    - b2 * (5.d0 + 2 * r2))) / (9 * x**3.d0)
                    
            gamma = bpr / (bmr * 3.d0 * x)
            gamma_r =  -o3 * (b * (3 * b2 - 5 * br) + r2 * r + b * (-2.d0 + r2)) &
                    / (bmr * bmr * ((1.d0 - bmr) * (1.d0 + bmr))**1.5)
            
            n = -4 * r * b / (bmr * bmr)
            n_r = -4 * b * bpr  / (bmr * bmr * bmr)
            
            m = 4 * r * b / ((1.d0 - bmr) * (1.d0 + bmr))
            m_r = (4 * b * ((1.d0 + b) * (1.d0 - b) + r2)) / ((bmr - 1.d0) * (bmr + 1.d0))**2.d0
            
            u = ((m - n) * (n - 1.d0) * ((m - 1.d0) * (m_r * alpha + 2 * m * alpha_r) - m_r * beta) &
                + m * (m_r - m_r * n + (m - 1.d0) * n_r) * gamma) &
                / (2 * (m - 1.d0) * m * (m - n) * (n - 1.d0))
            v = (2 * beta_r * m * (n - 1.d0) * n - m_r * (n - 1.d0) * n * (alpha + beta) + m * n_r * gamma) &
                / (2 * m * (n - 1.d0) * n)
            q = (2 * gamma_r * (m - n) * (n - 1.d0) * n - m * n_r * gamma &
                + n * (m_r - m_r * n + n * n_r) * gamma) / (2 * (m - n) * (n - 1.d0) * n)
            
            ome = Sqrt(1.d0 - m)
            o = 1.d0
            ellippi = cel((ome), 1.d0 - n, (o), (o))
            ellippi_r = ellippi * q
                
            eplusf = cel((ome), (o), alpha + beta, alpha * (1.d0 - m) + beta)
            eplusf_r = cel((ome), (o), u + v, u * (1.d0 - m) + v)
            
            Fl = eplusf + gamma * ellippi + pisixth * (1.d0 - Sign(1.d0, bmr))
            Fl_r = eplusf_r + ellippi_r
            goto 3
        end if
    end if
    
    if (b .eq. 0.d0) then
        Fl = pithird * (1.d0 - (1.d0 - r2) ** (1.5))
        Fl_r = -pi * r * Sqrt(1.d0 - r2)
        Fl_b = 0.d0
        goto 3
    end if
    
3   if (pflag) then
        Fcomplete(1) = (1.d0 - c1 - 2 * c2) * Fc + (c1 + 2 * c2) * Fl + c2 * Fq
        Fcomplete(2) = (1.d0 - c1 - 2 * c2) * Fc_b + (c1 + 2 * c2) * Fl_b + c2 * Fq_b
        Fcomplete(3) = (1.d0 - c1 - 2 * c2) * Fc_r + (c1 + 2 * c2) * Fl_r + c2 * Fq_r
        Fcomplete(4) = 0.d0
        Fcomplete(5) = 0.d0
        Fcomplete(6) = 0.d0
        Fcomplete(7) = - Fc + Fl 
        Fcomplete(8) = -2 * Fc + 2 * Fl + Fq
    else
        Fcomplete(1) = (1.d0 - c1 - 2 * c2) * Fc + (c1 + 2 * c2) * Fl + c2 * Fq
        Fcomplete(2) = 0.d0
        Fcomplete(3) = 0.d0
        Fcomplete(4) = (1.d0 - c1 - 2 * c2) * Fc_b + (c1 + 2 * c2) * Fl_b + c2 * Fq_b
        Fcomplete(5) = (1.d0 - c1 - 2 * c2) * Fc_r + (c1 + 2 * c2) * Fl_r + c2 * Fq_r
        Fcomplete(6) = 0.d0
        Fcomplete(7) = - Fc + Fl 
        Fcomplete(8) = -2 * Fc + 2 * Fl + Fq
    end if
    return

end function

function F(c1, c2, phi, r, b, phi_bp, phi_rp, phi_bm, phi_rm, phi_bpm, pflag)

    real*8, dimension(8) :: F
    
    logical :: pflag
    real*8 :: c1, c2, phi, cphi, r, b, Fc, Fq, Fl
    real*8 :: phi_bp, phi_rp, phi_bm, phi_rm, phi_bpm
    real*8 :: Fc_r, Fc_b, Fc_phi, Fq_r, Fq_b, Fq_phi, Fl_r, Fl_b, Fl_phi
    real*8 :: Fc_bp, Fc_rp, Fc_bm, Fc_rm, Fc_bpm
    real*8 :: Fq_bp, Fq_rp, Fq_bm, Fq_rm, Fq_bpm
    real*8 :: Fl_bp, Fl_rp, Fl_bm, Fl_rm, Fl_bpm
    real*8 :: gamma, d, s, n, m, x, y, ome, tans, sphi, br, bmr, bpr, o
    real*8 :: r2, b2, tanphi2, sinphi2, alpha, beta, alpha_r, beta_r, gamma_r
    real*8 :: d_phi, d_r, d_b, eplusf_phi, eplusf_r, eplusf_b, lpm, lpmo
    real*8 :: ellippi, eplusf, ellipf, eplusf_m
    real*8 :: m_r, m_b, n_r, n_b, ellippim, ellippi_m, ellippi_n, ellippi_r
    real*8 :: u, v, q, p
    
    r2 = r * r
    b2 = b * b
    bmr = b - r
    bpr = b + r
    br = b * r
    
    cphi = Cos(phi)
    sphi = Sin(phi)
    sinphi2 = Sin(phi * 0.5)
    tanphi2 = Tan(phi * 0.5)
    
    if (phi .eq. pi) then
        F = Fcomplete(c1, c2, r, b, pflag)
        return
    end if
        
    Fc = 0.5 * (r2 * phi - br * sphi)
    Fc_phi = 0.5 * (r2 - br * cphi)
    Fc_b = -0.5 * r * sphi
    Fc_r = 0.5 * (2 * r * phi - b * sphi)
    
    Fc_bp = Fc_phi * phi_bp
    Fc_rp = Fc_phi * phi_rp
    Fc_bm = Fc_phi * phi_bm
    Fc_rm = Fc_phi * phi_rm
    Fc_bpm = Fc_phi * phi_bpm 
    
    Fq = r2 * r2 * 0.25 * (phi - sphi * cphi * (4 * o3 * cphi * cphi - 1.d0)) &
        + b * (br * 0.5 * (phi * r - b * sphi * o3) &
        + r2 * r * sphi * 0.25 * ((4 * cphi * cphi - 1.d0) * o3 - 3.d0))

    Fq_phi = - o3 * 0.25 * (r * (b * (2 * b2 + 9 * r2) * cphi &
        - 3 * r * (2 * b2 + r2 + br * Cos(3 * phi)) + r2 * r * Cos(4 * phi)))

    Fq_r = b * (b * (r * phi - (b * sphi) * o3 * 0.5) &
        + 0.25 * r2 * (Sin(3 * phi) - 9 * sphi)) &
        + r2 * r* (phi - Sin(4*phi) * o3 * 0.25)

    Fq_b = b * (r2 * phi - br * sphi * 0.5) &
        + r2 * 0.25 * r * (Sin(3 * phi) * o3 - 3 * sphi)
                
    Fq_bp = Fq_phi * phi_bp
    Fq_rp = Fq_phi * phi_rp
    Fq_bm = Fq_phi * phi_bm
    Fq_rm = Fq_phi * phi_rm
    Fq_bpm = Fq_phi * phi_bpm 
        
    if (-c1 .eq. 2 * c2) then
        Fl = 0.d0
        Fl_phi = 0.d0
        Fl_r = 0.d0
        Fl_b = 0.d0
        goto 2
        
    else
        if (b == r) then
            if (r == 0.5) then
            
                Fl = phi * o3 * 0.5 + o3 * sinphi2 &
                    * (1.d0 - sinphi2 * sinphi2 * o3)
                Fl_phi = 0.5 * o3 * (1.d0 + 0.75 * Cos(phi * 0.5) + 0.25 * Cos(3 * phi * 0.5))
                Fl_r = 0.d0
                Fl_b = 0.d0
                goto 2
                
            else
                ! Still need to change this! 
                alpha = 4 * (2 * r2 - 1.d0) * o9
                alpha_r = 16 * r * o9
                
                beta = (1.d0 - 4 * r2) * o9
                beta_r = -8 * r * o9
                
                m = 4 * r2
                m_r = 8 * r
                
                d = phi * o3 * 0.5 - 2 * r2 * sphi * Sqrt(1.d0 + 2 * r2 * (cphi - 1.d0)) * o9
                d_phi = 0.5 * o3 - (r2 * ((2.d0 - 4 * r2) * cphi + r2 * (1.d0 + 3*Cos(2 * phi)))) &
                        / (9 * Sqrt(1.d0 - 2 * r2 + 2 * r2 * cphi))
                d_r = (-4 * r * (1.d0 - 3 * r2 + 3 * r2 * cphi) * sphi) &
                    / (9 * Sqrt(1.d0 - 2 * r2 + 2 * r2 * cphi))
                    
                Fl_phi = (alpha + beta - m * alpha * sphi * sphi) / Sqrt(1.d0 - m * sphi * sphi)
                
                eplusf = el2(Tan(phi * 0.5), Sqrt(1.d0 - m), alpha + beta, alpha * (1.d0 - m) + beta)
                eplusf_m = 0.5 * el2(Tan(phi * 0.5), Sqrt(1.d0 - m), &
                                     gamma - beta / m, gamma - alpha - beta / m)
                eplusf_r = el2(Tan(phi * 0.5), Sqrt(1.d0 - m), alpha_r + beta_r, alpha_r * (1.d0 - m) + beta_r)
                
                Fl = eplusf + d
                Fl_phi = Fl_phi + d_phi
                Fl_r = eplusf_m * m_r + eplusf_r + d_r
                goto 2
                
            end if
        else if (bpr == 1.d0) then
            ! also double check this at some point 
            y = Sqrt(br)
            Fl = phi * o3 * 0.5 + Atan((2 * r - 1.d0) / Tan(phi * 0.5)) * o3 &
                + Atan(2 * y * sinphi2 / (1.d0 - 2 * r)) * o3 &
                + pisixth *  Sqrt(1.d0 - 4 * r * b) / (2 * r - 1.d0) &
                + 2.d0 * o9 * y * (2 * r * (5 * r - 2.d0) &
                - 2 * br * cphi - 3.d0) * sinphi2
                
            Fl_r = ((3 * phi * (b + 2 * br - 1.d0)) / ((1.d0 - 2 * r)**2 * Sqrt(1 - 4 * br)) &
                + (phi / Tan(phi * 0.5)) / (1.d0 + (1 - 2*r)**2 / Tan(phi * 0.5)**2) &
                + (3 * b * (1.d0 + 2*r) * sinphi2) / (y * (1.d0 + 2 * r * (b + 2*r - 2.d0) - 2 * br * cphi)) &
                - (b * (3 + 2 * (6 - 25 * r) * r + 6 * br * cphi) * sphi) / y) * o9
                
            Fl_b = (r * ((3 * phi)/((1.d0 - 2 *r) * Sqrt(1 - 4 * br)) &
                + ((3*(1.d0 - 2*r)*sinphi2) / (1 + 2 * r * (b + 2*r - 2.d0) - 2 * br * cphi) &
                - (3 + 2 * (2 - 5*r) * r + 6 * br * cphi) * sphi)/y)) * o9
                
            Fl_phi = ((6 * Sqrt(1 - 4*br)) / (2*r - 1.d0) &
                + 2 * Atan((2*r - 1.d0) / Tan(phi * 0.5)) &
                - 8 * y * cphi * (3.d0 + 2 * (2.d0 - 5*r) * r + 2 * br * cphi) &
                + (12 * y *(2*r - 1.d0)*Cos(phi * 0.5)) / (-2 * r * (-2 + b + 2 * r - 2.d0) + 2 * br * cphi - 1.d0) &
                + (phi - 2 * phi * r) / (1.d0 + 2 * (r - 1.d0) * r &
                + 2 * (r - 1.d0) * r * cphi) + 16 * br**1.5 * sphi**2) * 0.25 * o9
            goto 2
            
        else if (bpr .gt. 1.d0) then
        
            y = Sqrt(br)
            
            alpha = 2 * y * (7 * r2 + b2 - 4.d0) * o9
            alpha_r = b * (b2 + 35 * r2 - 4.d0) / (9 * y)
            
            beta = 3.d0 + 2*r*(b2 * b + 5 * b2 * r + 3*r*(-2.d0 + r2) + b*(-4.d0 + 7*r2))
            beta = -beta * o9 * 0.5 / y
            beta_r = (b * (3.d0 - 2 * r * (b2 * b + 15 * b2 * r + 3 * r * (-6.d0 + 7 * r2) &
                + b * (-4.d0 + 35 * r2)))) / (36 * br**1.5)
            
            gamma = bpr * o3 / (2 * bmr * y)
            gamma_r = b * (4 * br + (r - b) * (r + b)) / (12 * bmr**2 * br**1.5)
            
            m = (1.d0 - bmr) * (1.d0 + bmr) / (4 * br)
            m_r = ((b + 1.d0) * (b - 1.d0) - r2) / (4 * br * r)
            
            n = 1.d0 - 1.d0 / (bmr * bmr)
            n_r = -2.d0 / (bmr**3.d0)
            
            u = ((m - n) * (n - 1.d0) * ((m - 1.d0) * (m_r * alpha + 2 * m * alpha_r) - m_r * beta) &
                + m * (m_r - m_r * n + (m - 1.d0) * n_r) * gamma) &
                / (2 * (m - 1.d0) * m * (m - n) * (n - 1.d0))
            v = (2 * beta_r * m * (n - 1.d0) * n - m_r * (n - 1.d0) * n * (alpha + beta) + m * n_r * gamma) &
                / (2 * m * (n - 1.d0) * n)
            q = (2 * gamma_r * (m - n) * (n - 1.d0) * n - m * n_r * gamma &
                + n * (m_r - m_r * n + n * n_r) * gamma) / (2 * (m - n) * (n - 1.d0) * n)
                
            ome = Sqrt(1.d0 - m)
            
            if (abs(sinphi2 / Sqrt(m) - 1.d0) .lt. 1.d-12) then
            
                d = phi * o3 * 0.5 - Atan2(bpr * sinphi2, bmr * Cos(phi * 0.5)) * o3
                d_phi = (r2 - br * cphi) / (3 * (b2 + r2 - 2 * br * cphi))
                d_r = - b * sphi / (3 * (b2 + r2 - 2 * br * cphi))
                
                Fl_phi = 0.d0
                p = 0.d0
                
                o = 1.d0
                ellippi = cel((ome), 1.d0 - n, (o), (o))
                ellippi_r = ellippi * q
                
                eplusf = cel((ome), (o), alpha + beta, alpha * (1.d0 - m) + beta)
                eplusf_r = cel((ome), (o), u + v, u * (1.d0 - m) + v)
                
            else
                ! Got work to do on this one because of the reciprocol modulus nonsense 
                
                y = Sqrt((1.d0 - b) * (1.d0 + b) - r2 + 2 * br * cphi)
                d = o3 * (phi * 0.5 - Atan2(bpr * sinphi2, bmr * Cos(phi * 0.5))) &
                    - (2 * br * o9) * sphi * y
                d_phi = o3 * o3 * r * (3 * (r - b * cphi) / (b2 + r2 - 2 * br * cphi) &
                    - 2 * b * cphi * y + 2 * b2 * r * sphi * sphi / y)
                d_r = o3 * o3 * b * (-3.d0 / (b2 + r2 - 2 * br * cphi) + 2 * (r2 - br * cphi) / y - 2 * y) * sphi
                
                Fl_phi = Cos(0.5 * phi) * (-(n * (alpha + 4 * beta)) &
                       + 4 * m * (alpha + 2*(beta + gamma)) &
                       + 4 * (m * alpha + n * beta) * cphi + n * alpha * Cos(2 * phi)) &
                       / (4 * Sqrt(m) * Sqrt(1.d0 + cphi) * Sqrt((-1.d0 + 2 * m + cphi) / m) * (2 * m - n + n * cphi))
                    
                p = - (4 * Sqrt(2.d0) * Cos(0.5 * phi)**2.d0 * (m_r * (-1.d0 + n) * (-2 * m + n) * (-m + n) &
                    * ((-1.d0 + m) * alpha - beta) + m * (-(m_r * (-1.d0 + n) * (-3 * n + 2 * m * (1.d0 + n))) &
                    + (-1.d0 + m) * (-1.d0 + 2 * m) * n * n_r) * gamma &
                    + n * (m_r * (m - n) * (-1.d0 + n) * ((-1.d0 + m)*alpha - beta) &
                    + m * (m_r - m_r * n + (-1.d0 + m) * n_r) * gamma) * cphi) * sinphi2) &
                    / (8 * (-1.d0 + m) * m**1.5 * (m - n) * (-1.d0 + n) * Sqrt(Cos(0.5 * phi)**2.d0) &
                    * Sqrt((-1.d0 + 2 * m + cphi) / m) * (2 * m - n + n * cphi))
                                        
                tans = 1.d0 / Sqrt(m / (sinphi2 * sinphi2) - 1.d0)
                
                ellippi = el3((tans), (ome), 1.d0 - n)
                ellippi_r = ellippi * q
                
                eplusf = el2((tans), (ome), alpha + beta, alpha * (1.d0 - m) + beta)
                eplusf_r = el2((tans), (ome), u + v, u * (1.d0 - m) + v)
                
            end if
            
            Fl = eplusf + gamma * ellippi + d
            Fl_phi = Fl_phi + d_phi  
            Fl_r = eplusf_r + ellippi_r + p + d_r
            goto 2
            
        else
        
            x = Sqrt((1.d0 - bmr) * (1.d0 + bmr))
            tans = tanphi2
            
            alpha = (7 * r2 + b2 - 4.d0) * x * o9
            alpha_r = (b * (b2 - 15 * br) + 3 * r * (6.d0 - 7 * r2) + b * (-4.d0 + 35 * r2)) &
                    / (9 * x)
            
            beta = (r2 * r2 + b2 * b2 + r2 - b2 * (5.d0 + 2 * r2) + 1.d0) / (9.d0 * x)
            beta_r = -o9 * (b2 * b2 * b - 5 * b2 * b2 * r + b2 * r * (11.d0 + 2 * r2) &
                   + b2 * b * (-5.d0 + 6 * r2) + b * (1.d0 - 3 * r2 - 7 * r2 * r2) + 3 * r * (-1.d0 - r2 + r2 * r2)) &
                   / ((1.d0 + b - r) * (1.d0 - b + r))**1.5
                    
            gamma = bpr / (bmr * 3.d0 * x)
            gamma_r =  -o3 * (b * (3 * b2 - 5 * br) + r2 * r + b * (-2.d0 + r2)) &
                    / (bmr * bmr * x**3.d0)
            
            n = -4 * r * b / (bmr * bmr)
            n_r = -4 * b * bpr / (bmr * bmr * bmr)
            
            m = 4 * r * b / ((1.d0 - bmr) * (1.d0 + bmr))
            m_r = (4 * b * ((1.d0 + b) * (1.d0 - b) + r2)) / ((bmr - 1.d0) * (bmr + 1.d0))**2.d0
            
            u = ((m - n) * (n - 1.d0) * ((m - 1.d0) * (m_r * alpha + 2 * m * alpha_r) - m_r * beta) &
                + m * (m_r - m_r * n + (m - 1.d0) * n_r) * gamma) &
                / (2 * (m - 1.d0) * m * (m - n) * (n - 1.d0))
            v = (2 * beta_r * m * (n - 1.d0) * n - m_r * (n - 1.d0) * n * (alpha + beta) + m * n_r * gamma) &
                / (2 * m * (n - 1.d0) * n)
            q = (2 * gamma_r * (m - n) * (n - 1.d0) * n - m * n_r * gamma &
                + n * (m_r - m_r * n + n * n_r) * gamma) / (2 * (m - n) * (n - 1.d0) * n)
            p = (-((m_r * (m - n) * (-2.d0 + n) * (-1.d0 + n) * beta + m * m_r * (-2.d0 + n) * (-1.d0 + n) * gamma &
                - (-2.d0 + m) * (-1.d0 + m) * n * n_r * gamma &
                + n * (-(m_r * (m - n) * (-1.d0 + n) * beta) &
                + m * (m_r - m_r * n + (-1.d0 + m) * n_r) * gamma) * cphi) * sphi)) &
                / (2 * (-1.d0 + m) * (m - n) * (-1.d0 + n) &
                * Sqrt(4.d0 - 2 * m + 2 * m * cphi) * (2.d0 - n + n * cphi))
            
            ome = Sqrt(1.d0 - m)
            
            y = Sqrt((1.d0 - b) * (1.d0 + b) - r2 + 2 * br * cphi)
            d = o3 * (phi * 0.5 - Atan((bpr / bmr) * tanphi2)) &
                - 2 * br * o9 * sphi * y
            d_r = o9 * b * sphi * (-3.d0 / (b2 + r2 - 2 * br * cphi) + 2 * (r2 - br * cphi) / y - 2 * y)
            d_phi = o9 * r * (3 * (r - b * cphi) / (b2 + r2 - 2 * br * cphi) &
                - 2 * b * cphi * y + 2 * b2 * r * sphi * sphi / y)
                
            ellippi = el3((tans), (ome), 1.d0 - n)
            ellippi_r = ellippi * q
                
            eplusf = el2((tans), (ome), alpha + beta, alpha * (1.d0 - m) + beta)
            eplusf_r = el2((tans), (ome), u + v, u * (1.d0 - m) + v)
            
            !Fl_phi = (-alpha - beta - gamma + (m * alpha + n * (alpha + beta)) * sphi * sphi &
            !        - m * n * alpha * sphi**4.d0) / (Sqrt(1.d0 - m * sphi * sphi) * (n * sphi * sphi - 1.d0))
            Fl_phi = (2 * alpha - m * alpha - n * alpha + 2 * beta - n * beta + 2 * gamma + m * alpha * cphi &
                    + n* alpha * cphi + n * beta * cphi + 2 * m * n * alpha * sinphi2**4) &
                    / (Sqrt(2.d0) * Sqrt(2.d0 - m + m * cphi) * (2.d0 - n + n * cphi))
            
            Fl = eplusf + gamma * ellippi + d
            Fl_phi = Fl_phi + d_phi
            Fl_r = eplusf_r + ellippi_r + p + d_r
            goto 2
            
        end if
        
        if (b .eq. 0.d0) then
            Fl = phi * (1.d0 - (1.d0 - r2) ** (1.5)) * o3
            Fl_phi = (1.d0 - (1.d0 - r2) ** (1.5)) * o3
            Fl_r = -pi * r * Sqrt(1.d0 - r2)
            goto 2
        end if
        
    end if
        
2   Fl_bp = Fl_phi * phi_bp
    Fl_rp = Fl_phi * phi_rp
    Fl_bm = Fl_phi * phi_bm
    Fl_rm = Fl_phi * phi_rm
    Fl_bpm = Fl_phi * phi_bpm 
    
    if (pflag) then
        Fq_bp = Fq_bp + Fq_b
        Fq_rp = Fq_rp + Fq_r
        Fc_bp = Fc_bp + Fc_b
        Fc_rp = Fc_rp + Fc_r
        Fl_bp = Fl_bp + Fl_b
        Fl_rp = Fl_rp + Fl_r
    else
        Fq_bm = Fq_bm + Fq_b
        Fq_rm = Fq_rm + Fq_r
        Fc_bm = Fc_bm + Fc_b
        Fc_rm = Fc_rm + Fc_r
        Fl_bm = Fl_bm + Fl_b
        Fl_rm = Fl_rm + Fl_r
    end if
    
    F(1) = (1.d0 - c1 - 2 * c2) * Fc + (c1 + 2 * c2) * Fl + c2 * Fq
    F(2) = (1.d0 - c1 - 2 * c2) * Fc_bp + (c1 + 2 * c2) * Fl_bp + c2 * Fq_bp
    F(3) = (1.d0 - c1 - 2 * c2) * Fc_rp + (c1 + 2 * c2) * Fl_rp + c2 * Fq_rp
    F(4) = (1.d0 - c1 - 2 * c2) * Fc_bm + (c1 + 2 * c2) * Fl_bm + c2 * Fq_bm
    F(5) = (1.d0 - c1 - 2 * c2) * Fc_rm + (c1 + 2 * c2) * Fl_rm + c2 * Fq_rm
    F(6) = (1.d0 - c1 - 2 * c2) * Fc_bpm + (c1 + 2 * c2) * Fl_bpm + c2 * Fq_bpm
    F(7) = - Fc + Fl 
    F(8) = -2 * Fc + 2 * Fl + Fq

    return
end function

end module phot