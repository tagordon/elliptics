module phot
use iso_c_binding
use ellip

implicit none

real*8, parameter :: pi = 3.14159265358979323846, pihalf = 1.5707963267948966, twopithree = 2.0943951023931953
real*8, parameter :: o3 = 0.33333333333333333333, o9 = 0.1111111111111111111, twopi = 6.283185307179586
real*8, parameter :: pithird = 1.0471975511965976, pisixth = 0.5235987755982988

contains

subroutine phis(rp, rm, bp, bm, bpm, theta, pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, &
                pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta)

    real*8 :: rp, rm, bp, bm, bpm, theta
    
    ! intersection angle from planet center, intersection from moon center, same 
    ! values relative to bp and bm vectors respectively
    real*8 :: pp, pm, pp1, pm1, pp2, pm2
    
    ! derivatives 
    real*8 :: pp_rp, pp_rm, pp_bpm ! pp_theta = 1
    real*8 :: pm_rp, pm_rm, pm_bpm, thetam_theta, thetam_bp, thetam_bpm
    
    ! Four times the area of the triangle formed by rm, rp, and bpm
    real*8 :: delta
    ! Variables used in sorting the sides of the triangle
    real*8 :: a, b, c, tmp
    
    ! angle between bpm vector and bm vector
    real*8 :: thetam
    
    a = bp
    b = bpm
    c = bm
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
    delta = Sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c)))
    thetam = Atan2(delta, (bm - bp) * (bm + bp) + bpm * bpm)
    
    ! find 4 * area of triangle using modified Heron's formula 
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
    delta = Sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c)))
    
    pm = Atan2(delta, (rm - rp) * (rm + rp) + bpm * bpm)   
    thetam_bp = bpm * Sin(theta) / (bm * bm)
    pm_bpm = ((rm + bpm) * (rm - bpm) - rp * rp) / (delta * bpm)
    pm_rp = 2 * rp / delta
    pm_rm = ((bpm - rm) * (bpm + rm) - rp * rp) / (delta * rm)
    thetam_theta = ((bpm - bm) * (bpm + bm) - bp * bp) / (2 * bm * bm)
    thetam_bpm = -bp * Sin(theta) / (bm * bm)
    
    pp = Atan2(delta, (rp - rm) * (rp + rm) + bpm * bpm)
    pp_bpm = ((rp - rm) * (rp + rm) - bpm * bpm) / (delta * bpm)
    pp_rp = ((bpm - rp) * (bpm + rp) - rm * rm) / (delta * rp)
    pp_rm = 2 * rm / delta
    
    ! this might be slower, but try it instead of the if-then statements  
    ! also check whether or not it's possible for pm1 or pp1 to be greater than pi
    pm1 = thetam + pm
    pm2 = thetam - pm
    pp1 = theta + pp
    pp2 = theta - pp
    
    !pm1 = pm1 - pi * (1.d0 - Sign(1.d0, pi - pm1))
    !pp1 = pp1 - pi * (1.d0 - Sign(1.d0, pi - pp1))
    
    if (pm1 .gt. pi) then
        pm1 = pm1 - twopi
    end if
    if (pm2 .gt. pi) then
        pm2 = pm2 - twopi
    end if
    if (pp1 .gt. pi) then
        pp1 = pp1 - twopi
    end if
    if (pp2 .gt. pi) then
        pp2 = pp2 - twopi
    end if

end

subroutine kappas_p(rp, bp, kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)

    ! kp = angle to intersection from center of planet, 
    ! kps = angle to intersection from center of star 
    real*8 :: rp, bp, kp, kps
    
    ! derivatives 
    real*8 :: kp_rp, kp_bp, kps_rp, kps_bp
    
    ! variables used in sorting sides of triangle
    real*8 :: a, b, c
    
    ! four times the area of the triangle with sides rp, bp, and 1
    real*8 :: delta
    
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
    delta = Sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c)))
    
    kps = Atan2(delta, (1.d0 - rp) * (1.d0 + rp) + bp * bp)
    kps_bp = ((1.d0 - bp) * (1.d0 + bp) - rp * rp) / (bp * delta)
    kps_rp = 2 * rp / delta
    
    kp = Atan2(delta, (rp - 1.d0) * (rp + 1.d0) + bp * bp)
    kp_bp = ((rp + bp) * (rp - bp) - 1.d0) / (bp * delta)
    kp_rp = ((bp + rp) * (bp - rp) - 1.d0) / (rp * delta)
end 

subroutine kappas_m(rm, bp, bm, bpm, theta, km, kms, km_rm, km_bp, km_bpm, km_theta, &
                    kms_rm, kms_bp, kms_bpm, kms_theta)

    ! km = angle to interection from center of moon, 
    ! kms = angle to intersection from center of planet
    real*8 :: rm, bp, bm, bpm, theta, km, kms
    
    ! derivatives
    real*8 :: km_rm, km_bp, km_bpm, km_theta, kms_rm, kms_bp, kms_bpm, kms_theta
    
    ! variables used in sorting sides of triangle
    real*8 :: a, b, c
    
    ! four times the area of the triangle with sides rm, bm, and 1
    real*8 :: delta
    
    ! some useful quantities
    real*8 :: denom, xs, xm, yp, ypm, ytheta
    
    xs = (1.d0 - bm) * (1.d0 + bm) - rm * rm
    xm = (rm - bm) * (rm + bm) - 1.d0
    yp = bp - bpm * Cos(theta)
    ypm = bpm - bp * Cos(theta)
    ytheta = bp * bpm * Sin(theta)
    
    if (bm .gt. 1.d0) then
        a = bm
        b = 1.d0
        c = rm
        if (rm .gt. bm) then
            b = rm
            c = bm
        end if
    else
        a = 1.d0
        b = bm
        c = rm
    end if
    delta = Sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c)))
    denom = 1.d0 / (delta * bm * bm)
    
    km = Atan2(delta, (rm - 1.d0) * (rm + 1.d0) + bm * bm)
    kms = Atan2(delta, (1.d0 - rm) * (1.d0 + rm) + bm * bm)
    
    km_rm = ((bm + rm) * (bm - rm) - 1.d0) / (rm * delta)
    kms_rm = 2 * rm * bm * bm * denom
    
    km_theta = ytheta * xm * denom
    kms_theta = ytheta * xs * denom
    
    km_bpm = ypm * xm * denom
    kms_bpm = ypm * xs * denom
    
    km_bp = yp * xm * denom
    kms_bp = yp * xs * denom
end 

subroutine bm_x(bp, bm, bpm, theta, bm_bp, bm_bpm, bm_theta)

    real*8 :: bp, bm, bpm, theta
    real*8 :: bm_bp, bm_bpm, bm_theta
    real*8 :: obm 
    
    obm = 1.d0 / bm
    bm_bp = (bp - bpm * Cos(theta)) * obm
    bm_bpm = (bpm - bp * Cos(theta)) * obm
    bm_theta = bp * bpm * Sin(theta) * obm
end 


! main loop to compute the flux at each timestep by finding the correct geometry and
! calling the integration routines 
subroutine flux(c1, c2, rp, rm, bp, bpm, theta, lc, j) bind(C, name="flux")

    integer (c_int), bind(C) :: j
    integer :: i
    real (c_double), bind(C) :: rp, rm
    real (c_double), bind(C), dimension(j) :: bp, theta, bpm
    real (c_double), bind(C), intent(out), dimension(8, j) :: lc
    real*8, dimension(8) :: f0
    real (c_double), bind(C) :: c1, c2
    real*8 :: of0
    
    ! half angles: from planet to planet/star intersection, moon to moon/star
    ! intersection, spanned by planet on limb of star, spanned by moon on 
    ! limb of star
    real*8 :: kp, km, kps, kms
        
    ! derivatives of above angles
    real*8 :: kp_rp, kp_bp, kps_rp, kps_bp
    real*8 :: km_rm, km_bp, km_bpm, km_theta
    real*8 :: kms_rm, kms_bp, kms_bpm, kms_theta
    
    ! angles to planet-moon intersection from planet center 
    !relative to bp vector (and derivatives)
    real*8 :: pp1, pp2
    real*8 :: pp_rp, pp_rm, pp_bpm ! pp_theta = 1
    
    ! angles to planet-moon intersection from moon center 
    ! relative to bm vector (and derivatves)
    real*8 :: pm1, pm2
    real*8 :: pm_rp, pm_rm, thetam_bp, pm_bpm, thetam_theta
    ! derivative of angle between bpm and bm vector with respect to bpm
    real*8 :: thetam_bpm
    
    ! used to determine cases for three body overlaps, might not be needed. 
    ! Check if some of these (costheta, cosphi) can be removed when optimizing things later 
    real*8 :: phi, phi_bpm, phi_bp, phi_bm, phi_theta, d1, d2, delta, a, b, c, tmp

    ! self explanatory 
    real*8 :: bpi, bmi, bpmi, rp2, rm2
    
    ! For chain rule stuff
    real*8, dimension(j) :: bm, ctheta
    real*8 :: obm, bm_bp, bm_bpm, bm_theta
    
    rp2 = rp * rp
    rm2 = rm * rm
    
    ctheta = Cos(theta)
    bm = Sqrt((bp - bpm)**2.d0 + 4 * bp * bpm * Sin(theta * 0.5)**2.d0)
    !obm = 1.d0 / bm
    !bm_bp = (bp - bpm * ctheta) * obms
    !bm_bpm = (bpm - bp * ctheta) * obm
    !bm_theta = bp * bpm * Sin(theta) * obm
    
    ! normalization factors 
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

        bpi = bp(i)
        bmi = bm(i)
        bpmi = bpm(i)
        
        if ((bpi .gt. rp + 1.d0) .AND. (bmi .gt. rm + 1.d0)) then
            ! neither planet nor moon overlap star
            lc(:, i) = f0 * of0
        else if (bpmi .gt. rp + rm) then
            if (bpi .gt. rp + 1.d0) then
                if (bmi .gt. rm + 1.d0) then
                    ! neither planet nor moon overlap star 
                    lc(:, i) = f0 * of0
                else
                    call bm_x(bpi, bmi, bpmi, theta(i), bm_bp, bm_bpm, bm_theta)
                    if (bmi + rm .le. 1.d0) then
                        ! moon completely overlaps star, planet is outside of star
                        lc(:, i) = (f0 - 2 * Fcomplete(c1, c2, rm, bmi, bm_bp, bm_bpm, bm_theta, .FALSE.)) * of0
                    else
                        ! moon partially overlaps star, planet is outside of star
                        call kappas_m(rm, bpi, bmi, bpmi, theta(i), km, kms, &
                                      km_rm, km_bp, km_bpm, km_theta, &
                                      kms_rm, kms_bp, kms_bpm, kms_theta)
                        lc(:, i) = 2 * (Fstar(c1, c2, pi - kms, 0.d0, -kms_rm, -kms_bp, -kms_bpm, -kms_theta) &
                                        - F(c1, c2, km, rm, bmi, 0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                            bm_bp, bm_bpm, bm_theta, .FALSE., .TRUE.)) * of0
                    end if
                end if
            else
                if (bmi .gt. rm + 1.d0) then
                    if (bpi + rp .le. 1.d0) then
                        ! planet completely overlaps star, moon is outside of star
                        lc(:, i) = (f0 - 2 * Fcomplete(c1, c2, rp, bpi, 0.d0, 0.d0, 0.d0, .TRUE.)) * of0
                    else
                        ! planet partially overlaps star, moon is outside of star
                        call kappas_p(rp, bpi, kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                        lc(:, i) = 2 * (Fstar(c1, c2, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) &
                                        - F(c1, c2, kp, rp, bpi, kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                            0.d0, 0.d0, 0.d0, .TRUE., .TRUE.)) * of0
                    end if
                else
                    call bm_x(bpi, bmi, bpmi, theta(i), bm_bp, bm_bpm, bm_theta)
                    if (bpi + rp .le. 1.d0) then
                        if (bmi + rm .le. 1.d0) then
                            ! moon and planet both completely overlap star, they do not overlap each othe
                            lc(:, i) = (f0 - 2 * (Fcomplete(c1, c2, rm, bmi, bm_bp, bm_bpm, bm_theta, .FALSE.) &
                                  + Fcomplete(c1, c2, rp, bpi, 0.d0, 0.d0, 0.d0, .TRUE.))) * of0
                        else
                            ! planet completely overlaps star, moon partially overlaps star, they do not overlap each other
                            call kappas_m(rm, bpi, bmi, bpmi, theta(i), km, kms, &
                                      km_rm, km_bp, km_bpm, km_theta, &
                                      kms_rm, kms_bp, kms_bpm, kms_theta)
                            lc(:, i) = 2 * (Fstar(c1, c2, pi - kms, 0.d0, -kms_rm, -kms_bp, -kms_bpm, -kms_theta) &
                                            - F(c1, c2, km, rm, bmi, 0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                                bm_bp, bm_bpm, bm_theta, .FALSE., .TRUE.) &
                                            - Fcomplete(c1, c2, rp, bpi, 0.d0, 0.d0, 0.d0, .TRUE.)) * of0
                        end if
                    else
                        if (bmi + rm .le. 1.d0) then
                            ! planet partially overlaps star, moon fully overlaps star, they do not overlap each other
                            call kappas_p(rp, bpi, kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                            lc(:, i) = 2 * (Fstar(c1, c2, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) & 
                                            - F(c1, c2, kp, rp, bpi, kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                0.d0, 0.d0, 0.d0, .TRUE., .TRUE.) &
                                            - Fcomplete(c1, c2, rm, bmi, bm_bp, bm_bpm, bm_theta, .FALSE.)) * of0
                        else
                            ! moon and planet both partially overlap star, but not each other
                            call kappas_p(rp, bpi, kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                            call kappas_m(rm, bpi, bmi, bpmi, theta(i), km, kms, &
                                      km_rm, km_bp, km_bpm, km_theta, &
                                      kms_rm, kms_bp, kms_bpm, kms_theta)
                            lc(:, i) = 2 * (Fstar(c1, c2, pi - (kps + kms), -kps_rp, -kms_rm, &
                                                  -(kps_bp + kms_bp), -kms_bpm, -kms_theta) &
                                            - F(c1, c2, kp, rp, bpi, kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                0.d0, 0.d0, 0.d0, .TRUE., .TRUE.) &
                                            - F(c1, c2, km, rm, bmi, 0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                                bm_bp, bm_bpm, bm_theta, .FALSE., .TRUE.)) * of0
                        end if
                    end if
                end if
            end if
        else
            if (bpi .gt. rp + 1.d0) then
                if (bmi .gt. rm + 1.d0) then
                    ! neither moon nor planet overlap star
                    lc(:, i) = f0
                else
                    ! moon partially overlaps star, planet does not overlap star
                    call bm_x(bpi, bmi, bpmi, theta(i), bm_bp, bm_bpm, bm_theta)
                    call kappas_m(rm, bpi, bmi, bpmi, theta(i), km, kms, &
                                      km_rm, km_bp, km_bpm, km_theta, &
                                      kms_rm, kms_bp, kms_bpm, kms_theta)
                    lc(:, i) = 2 * (Fstar(c1, c2, pi - kms, 0.d0, -kms_rm, -kms_bp, -kms_bpm, -kms_theta) &
                                    - F(c1, c2, km, rm, bmi, 0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                            bm_bp, bm_bpm, bm_theta, .FALSE., .TRUE.)) * of0

                end if
            else
                if (bmi .gt. rm + 1.d0) then
                    if (bpi + rp .le. 1.d0) then
                        ! planet fully overlaps star, moon does not overlap star
                        lc(:, i) = (f0 - 2 * Fcomplete(c1, c2, rp, bpi, 0.d0, 0.d0, 0.d0, .TRUE.)) * of0
                    else
                        ! planet partially overlaps star, moon is outside of star
                        call kappas_p(rp, bpi, kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                        lc(:, i) = 2 * (Fstar(c1, c2, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) &
                                        - F(c1, c2, kp, rp, bpi, kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                            0.d0, 0.d0, 0.d0, .TRUE., .TRUE.)) * of0
                    end if
                else
                    if (bpi + rp .le. 1.d0) then
                        if (bmi + rm .le. 1.d0) then
                            if (bpmi + rm .le. rp) then
                                ! moon and planet both overlap star, moon fully overlapped by planet
                                lc(:, i) = (f0 - 2 * Fcomplete(c1, c2, rp, bpi, 0.d0, 0.d0, 0.d0, .TRUE.)) * of0
                            else
                                ! moon and planet both overlap star, moon and planet partially overlap each other 
                                ! bookmark
                                call bm_x(bpi, bmi, bpmi, theta(i), bm_bp, bm_bpm, bm_theta)
                                call phis(rp, rm, bpi, bmi, bpmi, theta(i), pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, &
                                          pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta)
                                lc(:, i) = (f0 - Arc(c1, c2, pp1, pp2, rp, bpi, &
                                                     pp_rp, pp_rm, 0.d0, pp_bpm, 1.d0, &
                                                     -pp_rp, -pp_rm, 0.d0, -pp_bpm, 1.d0, &
                                                     0.d0, 0.d0, 0.d0, .TRUE., .FALSE., .FALSE.) &
                                                - Arc(c1, c2, pm1, pm2, rm, bmi, &
                                                     pm_rp, pm_rm, thetam_bp, pm_bpm + thetam_bpm, thetam_theta, &
                                                     -pm_rp, -pm_rm, thetam_bp, -pm_bpm + thetam_bpm, thetam_theta, &
                                                     bm_bp, bm_bpm, bm_theta, .FALSE., .FALSE., .FALSE.)) * of0
                            end if
                        else
                            call bm_x(bpi, bmi, bpmi, theta(i), bm_bp, bm_bpm, bm_theta)
                            call phis(rp, rm, bpi, bmi, bpmi, theta(i), pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, &
                                      pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta)
                            call kappas_m(rm, bpi, bmi, bpmi, theta(i), km, kms, &
                                      km_rm, km_bp, km_bpm, km_theta, &
                                      kms_rm, kms_bp, kms_bpm, kms_theta)
                            lc(:, i) = (2 * Fstar(c1, c2, pi - kms, 0.d0, -kms_rm, -kms_bp, -kms_bpm, -kms_theta) &
                                        - Arc(c1, c2, -km, pm2, rm, bmi, &
                                              0.d0, -km_rm, -km_bp, -km_bpm, -km_theta, &
                                              -pm_rp, -pm_rm, thetam_bp, -pm_bpm + thetam_bpm, thetam_theta, &
                                              bm_bp, bm_bpm, bm_theta, .FALSE., .TRUE., .FALSE.) &
                                        - Arc(c1, c2, pm1, km, rm, bmi, &
                                              pm_rp, pm_rm, thetam_bp, pm_bpm + thetam_bpm, thetam_theta, &
                                              0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                              bm_bp, bm_bpm, bm_theta, .FALSE., .FALSE., .TRUE.) &
                                        - Arc(c1, c2, pp1, pp2, rp, bpi, &
                                              pp_rp, pp_rm, 0.d0, pp_bpm, 1.d0, &
                                              -pp_rp, -pp_rm, 0.d0, -pp_bpm, 1.d0, &
                                              0.d0, 0.d0, 0.d0, .TRUE., .FALSE., .FALSE.)) * of0
                        end if
                    else
                        if (bmi + rm .le. 1.d0) then 
                            if (bpmi + rm .le. rp) then
                                ! planet partially overlaps star, moon fully overlaps star but is completely overlapped by planet 
                                call kappas_p(rp, bpi, kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                                lc(:, i) = 2 * (Fstar(c1, c2, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) &
                                                - F(c1, c2, kp, rp, bpi, kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                    0.d0, 0.d0, 0.d0, .TRUE., .TRUE.)) * of0
                            else
                                ! planet partially overlaps star, moon fully overlaps star and only partially overlaps planet
                                call bm_x(bpi, bmi, bpmi, theta(i), bm_bp, bm_bpm, bm_theta)
                                call phis(rp, rm, bpi, bmi, bpmi, theta(i), pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, &
                                      pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta)
                                call kappas_p(rp, bpi, kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                                lc(:, i) = (2 * Fstar(c1, c2, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) &
                                            - Arc(c1, c2, -kp, pp2, rp, bpi, &
                                                  -kp_rp, 0.d0, -kp_bp, 0.d0, 0.d0, &
                                                  -pp_rp, -pp_rm, 0.d0, -pp_bpm, 1.d0, &
                                                  0.d0, 0.d0, 0.d0, .TRUE., .TRUE., .FALSE.) &
                                            - Arc(c1, c2, pp1, kp, rp, bpi, &
                                                  pp_rp, pp_rm, 0.d0, pp_bpm, 1.d0, &
                                                  kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                  0.d0, 0.d0, 0.d0, .TRUE., .FALSE., .TRUE.) &
                                            - Arc(c1, c2, pm1, pm2, rm, bmi, &
                                                  pm_rp, pm_rm, thetam_bp, pm_bpm + thetam_bpm, thetam_theta, &
                                                  -pm_rp, -pm_rm, thetam_bp, -pm_bpm + thetam_bpm, thetam_theta, &
                                                  bm_bp, bm_bpm, bm_theta, .FALSE., .FALSE., .FALSE.)) * of0
                            end if
                        else
                            if (bpmi + rm .le. rp) then
                                ! planet and moon both partially overlap star but moon is fully overlapped by the planet
                                call kappas_p(rp, bpi, kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                                lc(:, i) = 2 * (Fstar(c1, c2, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) &
                                                - F(c1, c2, kp, rp, bpi, kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                    0.d0, 0.d0, 0.d0, .TRUE., .TRUE.)) * of0
                            else
                                !call compute_theta(rp,  bpi, theta, phip, theta_bp, theta_rp, phip_bp, phip_rp)
                                !call compute_theta(rm,  bmi, theta, phim, theta_bm, theta_rm, phim_bm, phim_rm)
                                call kappas_p(rp, bpi, kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                                call kappas_m(rm, bpi, bmi, bpmi, theta(i), km, kms, &
                                      km_rm, km_bp, km_bpm, km_theta, &
                                      kms_rm, kms_bp, kms_bpm, kms_theta)
                                
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
                                
                                delta = Sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c)))
                                phi = Atan2(delta, (bmi - bpmi) * (bmi + bpmi) + bpi * bpi)
                                
                                ! Probably need phi_theta rather than phi_bm
                                phi_bpm = bpi * Sin(theta(i)) / (bmi * bmi)
                                phi_theta = bpmi * (bpi * Cos(theta(i)) - bpmi) / (bmi * bmi)
                                phi_bp = - bpmi * Sin(theta(i)) / (bmi * bmi)
                                
                                call phis(rp, rm, bpi, bmi, bpmi, theta(i), pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, &
                                          pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta)
                                
                                if (phi + kms .le. kps) then         
                                        if (pp2 .gt. kp) then
                                            ! planet and moon both partially overlap the star and each other but the 
                                            ! moon-star overlap is contained within the planet-star overlap
                                            lc(:, i) = 2 * (Fstar(c1, c2, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) &
                                                            - F(c1, c2, kp, rp, bpi, kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                                0.d0, 0.d0, 0.d0, .TRUE., .TRUE.)) * of0
                                        else
                                            ! planet and moon both partially overlap star and each other but the 
                                            ! planet-star intersections are overlapped by the planet
                                            call bm_x(bpi, bmi, bpmi, theta(i), bm_bp, bm_bpm, bm_theta)
                                            lc(:, i) = (2 * Fstar(c1, c2, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) &
                                                        - Arc(c1, c2, -kp, pp2, rp, bpi, &
                                                              -kp_rp, 0.d0, -kp_bp, 0.d0, 0.d0, &
                                                              -pp_rp, -pp_rm, 0.d0, -pp_bpm, 1.d0, &
                                                              0.d0, 0.d0, 0.d0, .TRUE., .TRUE., .FALSE.) &
                                                        - Arc(c1, c2, pp1, kp, rp, bpi, &
                                                              pp_rp, pp_rm, 0.d0, pp_bpm, 1.d0, &
                                                              kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                              0.d0, 0.d0, 0.d0, .TRUE., .FALSE., .TRUE.) &
                                                        - Arc(c1, c2, pm1, pm2, rm, bmi, &
                                                              pm_rp, pm_rm, thetam_bp, pm_bpm + thetam_bpm, thetam_theta, &
                                                              -pm_rp, -pm_rm, thetam_bp, -pm_bpm + thetam_bpm, thetam_theta, &
                                                              bm_bp, bm_bpm, bm_theta, .FALSE., .FALSE., .FALSE.)) * of0
                                        end if
                                else if (phi + kps .le. kms) then
                                    call bm_x(bpi, bmi, bpmi, theta(i), bm_bp, bm_bpm, bm_theta)
                                    if ((bpi - rp) .le. (bmi - rm)) then
                                        ! planet and moon both partially overlap the star and each other but the 
                                        ! planet-star intersections are overlapped by the moon
                                        ! I'm not sure this is physical either -- can you draw a diagram where 
                                        ! the moon overlaps both of the planet-star 
                                        ! intersections without the planet-star overlap being entirely within 
                                        ! the moon-star region of overlap?
                                        lc(:, i) = (2 * Fstar(c1, c2, pi - kms, 0.d0, -kms_rm, -kms_bp, -kms_bpm, -kms_theta) &
                                                        - Arc(c1, c2, -km, pm2, rm, bmi, &
                                                              0.d0, -km_rm, -km_bp, -km_bpm, -km_theta, &
                                                              -pm_rp, -pm_rm, thetam_bp, -pm_bpm + thetam_bpm, thetam_theta, &
                                                              bm_bp, bm_bpm, bm_theta, .FALSE., .TRUE., .FALSE.) &
                                                        - Arc(c1, c2, pm1, km, rm, bmi, &
                                                              pm_rp, pm_rm, thetam_bp, pm_bpm + thetam_bpm, thetam_theta, &
                                                              km_rm, 0.d0, km_bp, km_bpm, km_theta, &
                                                              bm_bp, bm_bpm, bm_theta, .FALSE., .FALSE., .TRUE.) &
                                                        - Arc(c1, c2, pp1, pp2, rp, bpi, &
                                                              pp_rp, pp_rm, 0.d0, pp_bpm, 1.d0, &
                                                              -pp_rp, -pp_rm, 0.d0, -pp_bpm, 1.d0, &
                                                              0.d0, 0.d0, 0.d0, .TRUE., .FALSE., .FALSE.)) * of0
                                    else
                                        ! planet and moon both partially overlap the star and each other but 
                                        ! the planet-star overlap is  entirely within the moon-star overlap
                                        lc(:, i) = 2 * (Fstar(c1, c2, pi - kms, 0.d0, -kms_rm, -kms_bp, -kms_bpm, -kms_theta) &
                                                        - F(c1, c2, km, rm, bmi, 0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                                            bm_bp, bm_bpm, bm_theta, .FALSE., .TRUE.)) * of0
                                    end if
                                else
                                    ! bookmark
                                    d1 = rm2 + bmi * bmi - 2 * rm * bmi * Cos(pm2)
                                    d2 = rm2 + bmi * bmi - 2 * rm * bmi * Cos(pm1)
                                    call bm_x(bpi, bmi, bpmi, theta(i), bm_bp, bm_bpm, bm_theta)
                                    if (d1 .gt. 1.d0) then
                                        ! planet and moon both partially overlap star and each other, 
                                        ! but the planet/moon overlap does not overlap the star
                                        lc(:, i) = 2 * (Fstar(c1, c2, pi - (kps + kms), -kps_rp, -kms_rm, &
                                                              -(kps_bp + kms_bp), -kms_bpm, -kms_theta) &
                                                        - F(c1, c2, kp, rp, bpi, kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                            0.d0, 0.d0, 0.d0, .TRUE., .TRUE.) &
                                                        - F(c1, c2, km, rm, bmi, 0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                                            bm_bp, bm_bpm, bm_theta, .FALSE., .TRUE.)) * of0
                                    else if (d2 .le. 1.d0) then
                                        ! planet and moon both partially overlap star and each other, 
                                        ! with the planet/moon overlap fully overlapping the star
                                        lc(:, i) = (2 * Fstar(c1, c2, pi - (kps + kms), -kps_rp, -kms_rm, &
                                                              -(kps_bp + kms_bp), -kms_bpm, -kms_theta) &
                                                        - Arc(c1, c2, -km, -pm1, rm, bmi, &
                                                              0.d0, -km_rm, -km_bp, -km_bpm, -km_theta, &
                                                              -pm_rp, -pm_rm, -thetam_bp, -pm_bpm - thetam_bpm, -thetam_theta, &
                                                              bm_bp, bm_bpm, bm_theta, .FALSE., .TRUE., .FALSE.) &
                                                        - Arc(c1, c2, -pm2, km, rm, bmi, &
                                                              pm_rp, pm_rm, -thetam_bp, pm_bpm - thetam_bpm, -thetam_theta, &
                                                              0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                                              bm_bp, bm_bpm, bm_theta, .FALSE., .FALSE., .TRUE.) &
                                                        - Arc(c1, c2, pp1, kp, rp, bpi, &
                                                              pp_rp, pp_rm, 0.d0, pp_bpm, 1.d0, &
                                                              kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                              0.d0, 0.d0, 0.d0, .TRUE., .FALSE., .TRUE.) &
                                                        - Arc(c1, c2, -kp, pp2, rp, bpi, &
                                                              -kp_rp, 0.d0, -kp_bp, 0.d0, 0.d0, &
                                                              -pp_rp, -pp_rm, 0.d0, -pp_bpm, 1.d0, &
                                                              0.d0, 0.d0, 0.d0, .TRUE., .TRUE., .FALSE.)) * of0
                                    else
                                        ! planet and moon both partially overlap star and each other, 
                                        ! with the planet/moon overlap partially overlapping the star
                                        lc(:, i) = (2 * Fstar(c1, c2, pi - 0.5 * (kps + kms + phi), -0.5 * kps_rp, -0.5 * kms_rm, &
                                                              -0.5 * (kps_bp + kms_bp + phi_bp), -0.5 * (kms_bpm + phi_bpm), &
                                                              -0.5 * (kms_theta + phi_theta)) &
                                                        - Arc(c1, c2, -pm2, km, rm, bmi, &
                                                              pm_rp, pm_rm, -thetam_bp, pm_bpm - thetam_bpm, -thetam_theta, &
                                                              0.d0, km_rm, km_bp, km_bpm, km_theta,  &
                                                              bm_bp, bm_bpm, bm_theta, .FALSE., .FALSE., .TRUE.) &
                                                        - Arc(c1, c2, -kp, pp2, rp, bpi, &
                                                              -kp_rp, 0.d0, -kp_bp, 0.d0, 0.d0, &
                                                              -pp_rp, -pp_rm, 0.d0, -pp_bpm, 1.d0, &
                                                              0.d0, 0.d0, 0.d0, .TRUE., .TRUE., .FALSE.)) * of0
                                    end if
                                end if
                            end if
                        end if
                    end if
                end if
            end if  
        end if
        lc(:, i) = lc(:, i) - f0 * of0
1   end do
    return
    
end

! work out the right sign and order of the integration and call the integration routine 
! to integrate along an arbitrary arc of the planet or moon 
function Arc(c1, c2, phi1, phi2, r, b, phi1_rp, phi1_rm, phi1_bp, phi1_bpm, phi1_theta, &
            phi2_rp, phi2_rm, phi2_bp, phi2_bpm, phi2_theta, bm_bp, bm_bpm, bm_theta, pflag, limbflag1, limbflag2)
                    
    real*8, dimension(8) :: Arc

    logical :: pflag, limbflag1, limbflag2
    real*8 :: phi1, phi2, r, b, c1, c2, bm_bp, bm_bpm, bm_theta
    real*8 :: phi1_rp, phi1_rm, phi1_bp, phi1_bpm, phi1_theta
    real*8 :: phi2_rp, phi2_rm, phi2_bp, phi2_bpm, phi2_theta
    real*8 :: const, lin, quad
        
    if (phi1 < 0) then
        if (phi2 > 0) then
            Arc = F(c1, c2, phi2, r, b, phi2_rp, phi2_rm, phi2_bp, phi2_bpm, phi2_theta, &
                    bm_bp, bm_bpm, bm_theta, pflag, limbflag2) &
                + F(c1, c2, -phi1, r, b, -phi1_rp, -phi1_rm, -phi1_bp, -phi1_bpm, -phi1_theta, &
                    bm_bp, bm_bpm, bm_theta, pflag, limbflag1)
            return
        else
            if (phi2 < phi1) then
                Arc = 2 * Fcomplete(c1, c2, r, b, bm_bp, bm_bpm, bm_theta, pflag) &
                    + F(c1, c2, -phi1, r, b, -phi1_rp, -phi1_rm, -phi1_bp, -phi1_bpm, -phi1_theta, &
                        bm_bp, bm_bpm, bm_theta, pflag, limbflag1) &
                    - F(c1, c2, -phi2, r, b, -phi2_rp, -phi2_rm, -phi2_bp, -phi2_bpm, -phi2_theta, &
                        bm_bp, bm_bpm, bm_theta, pflag, limbflag2)
                return
            else
                Arc = - F(c1, c2, -phi2, r, b, -phi2_rp, -phi2_rm, -phi2_bp, -phi2_bpm, -phi2_theta, &
                          bm_bp, bm_bpm, bm_theta, pflag, limbflag2) &
                      + F(c1, c2, -phi1, r, b, -phi1_rp, -phi1_rm, -phi1_bp, -phi1_bpm, -phi1_theta, &
                          bm_bp, bm_bpm, bm_theta, pflag, limbflag1)
                return
            end if
        end if
    else
        if (phi2 < 0) then
            Arc = 2 * Fcomplete(c1, c2, r, b, bm_bp, bm_bpm, bm_theta, pflag) &
                - F(c1, c2, phi1, r, b, phi1_rp, phi1_rm, phi1_bp, phi1_bpm, phi1_theta, &
                    bm_bp, bm_bpm, bm_theta, pflag, limbflag1) &
                - F(c1, c2, -phi2, r, b, -phi2_rp, -phi2_rm, -phi2_bp, -phi2_bpm, -phi2_theta, &
                    bm_bp, bm_bpm, bm_theta, pflag, limbflag2)
            return
        else
            if (phi2 < phi1) then
                Arc = 2 * Fcomplete(c1, c2, r, b, bm_bp, bm_bpm, bm_theta, pflag) &
                    + F(c1, c2, phi2, r, b, phi2_rp, phi2_rm, phi2_bp, phi2_bpm, phi2_theta, &
                        bm_bp, bm_bpm, bm_theta, pflag, limbflag2) &
                    - F(c1, c2, phi1, r, b, phi1_rp, phi1_rm, phi1_bp, phi1_bpm, phi1_theta, &
                        bm_bp, bm_bpm, bm_theta, pflag, limbflag1)
                return
            else
                Arc = F(c1, c2, phi2, r, b, phi2_rp, phi2_rm, phi2_bp, phi2_bpm, phi2_theta, &
                        bm_bp, bm_bpm, bm_theta, pflag, limbflag2) &
                    - F(c1, c2, phi1, r, b, phi1_rp, phi1_rm, phi1_bp, phi1_bpm, phi1_theta, &
                        bm_bp, bm_bpm, bm_theta, pflag, limbflag1)
                return
            end if
        end if
    end if
    
    return
    
end function

! integrate along the limb of the star
function Fstar(c1, c2, phi, phi_rp, phi_rm, phi_bp, phi_bpm, phi_theta)

    real*8, dimension(8) :: Fstar

    real*8 :: c1, c2, phi, Fc, Fq, Fl
    real*8 :: Fc_bp, Fc_rp, Fc_bm, Fc_rm, Fc_bpm, Fc_theta
    real*8 :: Fq_bp, Fq_rp, Fq_bm, Fq_rm, Fq_bpm, Fq_theta
    real*8 :: Fl_bp, Fl_rp, Fl_bm, Fl_rm, Fl_bpm, Fl_theta
    real*8 :: Fc_phi, Fq_phi, Fl_phi
    real*8 :: phi_bp, phi_rp, phi_bm, phi_rm, phi_bpm, phi_theta
    real*8 :: cc, cl, cq
    
    Fc = 0.5 * phi
    Fc_phi = 0.5
    Fc_rp = Fc_phi * phi_rp
    Fc_rm = Fc_phi * phi_rm
    Fc_bp = Fc_phi * phi_bp
    Fc_bpm = Fc_phi * phi_bpm
    Fc_theta = Fc_phi * phi_theta

    Fq = 0.25 * (phi + 0.5 * o3 * (Sin(2 * phi) - Sin(4 * phi)))
    Fq_phi = 0.5 * o3 * Cos(phi)**2.d0 * (5.d0 - 4 * Cos(2 * phi))
    Fq_rp = Fq_phi * phi_rp
    Fq_rm = Fq_phi * phi_rm
    Fq_bp = Fq_phi * phi_bp
    Fq_bpm = Fq_phi * phi_bpm
    Fq_theta = Fq_phi * phi_theta
    
    Fl = phi * o3
    Fl_phi = o3
    Fl_rp = Fl_phi * phi_rp
    Fl_rm = Fl_phi * phi_rm
    Fl_bp = Fl_phi * phi_bp
    Fl_bpm = Fl_phi * phi_bpm
    Fl_theta = Fl_phi * phi_theta
    
    cc = 1.d0 - c1 - 2 * c2
    cl = c1 + 2 * c2
    cq = c2
        
    Fstar(1) = cc * Fc + cl * Fl + cq * Fq
    Fstar(2) = cc * Fc_rp + cl * Fl_rp + cq * Fq_rp
    Fstar(3) = cc * Fc_rm + cl * Fl_rm + cq * Fq_rm
    Fstar(4) = cc * Fc_bp + cl * Fl_bp + cq * Fq_bp
    Fstar(5) = cc * Fc_bpm + cl * Fl_bpm + cq * Fq_bpm
    Fstar(6) = cc * Fc_theta + cl * Fl_theta + cq * Fq_theta
    Fstar(7) = - Fc + Fl 
    Fstar(8) = -2 * Fc + 2 * Fl + Fq
    return
    
end function

! integrate around the entire planet/moon 
function Fcomplete(c1, c2, r, b, bm_bp, bm_bpm, bm_theta, pflag)

    real*8, dimension(8) :: Fcomplete
    
    ! Are we integrating along the edge of the moon or the planet? 
    logical :: pflag
    
    ! Limb darkening params
    real*8 :: c1, c2, cc, cl, cq
    
    ! self explanatory
    real*8 :: r, b
    
    ! derivatives of input parameters 
    real*8 :: bm_bp, bm_bpm, bm_theta
    
    ! convenient parameters
    real*8 :: r2, b2, br, bmr, bpr, obmr
    real*8 :: x, y, ox, ome, o
    
    ! components of flux and their derivatives
    real*8 :: Fc, Fc_r, Fc_b, Fc_bp, Fc_bpm, Fc_theta
    real*8 :: Fq, Fq_r, Fq_b, Fq_bp, Fq_bpm, Fq_theta
    real*8 :: Fl, Fl_r, Fl_b, Fl_bp, Fl_bpm, Fl_theta
    
    ! For the integral
    real*8 :: alpha, beta, gamma, d, n, m
    real*8 :: d_r, d_b, sgn
    real*8 :: ur, vr, ub, vb, sqomm
    
    ! Elliptic integrals
    real*8 :: ellippi, ellipe, ellipf
    real*8 :: eplusf, eplusf_r, eplusf_b
    
    r2 = r * r
    b2 = b * b
    bmr = b - r
    obmr = 1.d0 / bmr
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
    else
        x = Sqrt((1.d0 - bmr) * (1.d0 + bmr))
        ox = 1.d0 / x
            
        alpha = (7 * r2 + b2 - 4.d0) * o9 * x
        beta = (r2 * r2 + b2 * b2 + r2 - b2 * (5.d0 + 2 * r2) + 1.d0) * o9 * ox        
            
        n = -4 * br * obmr * obmr
        m = 4 * br * ox * ox
            
        ur = 2 * r * x
        ub = x * ((b + 1.d0) * (b - 1.d0) + r2) / (3 * b)
        vb = ((-1.d0 + b2)**2.d0 - 2 * (1.d0 + b2) * r2 + r2 * r2) * o3 * ox / b
            
        sqomm = Sqrt(1.d0 - m)
        o = 1.d0
            
        if (b .eq. r) then 
            ellippi = 0.d0
            sgn = 0.d0
            gamma = 0.d0
        else
            ellippi = cel((sqomm), (1.d0 - n), (o), (o))
            sgn = Sign(1.d0, bmr)
            gamma = bpr * ox * o3 * obmr
        end if 
            
        ellipe = cel((sqomm), (o), (o), 1.d0 - m)
        ellipf = cel((sqomm), (o), (o), (o))
             
        eplusf = alpha * ellipe + beta * ellipf
        eplusf_r = ur * ellipe
        eplusf_b = ub * ellipe + vb * ellipf
            
        Fl = eplusf + gamma * ellippi + pisixth * (1.d0 - sgn)
        Fl_r = eplusf_r
        Fl_b = eplusf_b
    end if
    
    if (b .eq. 0.d0) then
        Fl = pithird * (1.d0 - (1.d0 - r2) ** (1.5))
        Fl_r = -pi * r * Sqrt(1.d0 - r2)
        Fl_b = 0.d0
    end if

    cc = 1.d0 - c1 - 2 * c2
    cl = c1 + 2 * c2
    cq = c2
    Fcomplete = 0.d0
    
    if (pflag) then
        Fcomplete(1) = cc * Fc + cl * Fl + cq * Fq
        Fcomplete(2) = cc * Fc_r + cl * Fl_r + cq * Fq_r
        Fcomplete(4) = cc * Fc_b + cl * Fl_b + cq * Fq_b
        Fcomplete(7) = - Fc + Fl 
        Fcomplete(8) = -2 * Fc + 2 * Fl + Fq
    else
        Fcomplete(1) = cc * Fc + cl * Fl + cq * Fq
        Fcomplete(3) = cc * Fc_r + cl * Fl_r + cq * Fq_r
        Fcomplete(4) = (cc * Fc_b + cl * Fl_b + cq * Fq_b) * bm_bp
        Fcomplete(5) = (cc * Fc_b + cl * Fl_b + cq * Fq_b) * bm_bpm
        Fcomplete(6) = (cc * Fc_b + cl * Fl_b + cq * Fq_b) * bm_theta
        Fcomplete(7) = - Fc + Fl 
        Fcomplete(8) = -2 * Fc + 2 * Fl + Fq
    end if
    return

end function

! evaluate the integral at one arbitrary limit along the planet or moon's boundary 
function F(c1, c2, phi, r, b, phi_rp, phi_rm, phi_bp, phi_bpm, phi_theta, bm_bp, bm_bpm, bm_theta, pflag, limbflag)

    real*8, dimension(8) :: F
    
    ! Are we integrating along the edge of the moon or the planet? 
    logical :: pflag, limbflag
    
    ! Limb darkening params
    real*8 :: c1, c2, cc, cl, cq
    
    ! self explanatory
    real*8 :: phi, r, b
    
    ! derivatives of input parameters 
    real*8 :: phi_bp, phi_rp, phi_bm, phi_rm, phi_bpm, phi_theta, bm_bp, bm_bpm, bm_theta
    
    ! convenient parameters
    real*8 :: sphi, cphi, tphihalf, sphihalf, cphihalf
    real*8 :: r2, b2, br, bmr, bpr, obmr
    real*8 :: x, y, z, ox, oy, oz, ome, tans, o
    
    ! components of flux and their derivatives
    real*8 :: Fc, Fc_phi, Fc_r, Fc_b, Fc_rp, Fc_rm, Fc_bp, Fc_bm, Fc_bpm, Fc_theta
    real*8 :: Fq, Fq_phi, Fq_r, Fq_b, Fq_rp, Fq_rm, Fq_bp, Fq_bm, Fq_bpm, Fq_theta
    real*8 :: Fl, Fl_phi, Fl_r, Fl_b, Fl_rp, Fl_rm, Fl_bp, Fl_bm, Fl_bpm, Fl_theta
    
    ! For the integral
    real*8 :: alpha, beta, gamma, d, n, m
    real*8 :: d_phi, d_r, d_b
    real*8 :: ur, vr, ub, vb, pr, pb, sqomm
    
    ! Elliptic integrals
    real*8 :: ellippi, ellipe, ellipf
    real*8 :: eplusf, eplusf_r, eplusf_b
    
    r2 = r * r
    b2 = b * b
    bmr = b - r
    obmr = 1.d0 / bmr
    bpr = b + r
    br = b * r
    
    cphi = Cos(phi)
    sphi = Sin(phi)
    sphihalf = Sin(phi * 0.5)
    cphihalf = Cos(phi * 0.5)
    tphihalf = sphihalf / cphihalf
    
    if (phi .eq. pi) then
        F = Fcomplete(c1, c2, r, b, bm_bp, bm_bpm, bm_theta, pflag)
        return
    end if
        
    Fc = 0.5 * (r2 * phi - br * sphi)
    Fc_phi = 0.5 * (r2 - br * cphi)
    Fc_b = -0.5 * r * sphi
    Fc_r = 0.5 * (2 * r * phi - b * sphi)
    
    Fc_rp = Fc_phi * phi_rp
    Fc_rm = Fc_phi * phi_rm
    Fc_bp = Fc_phi * phi_bp
    Fc_bpm = Fc_phi * phi_bpm
    Fc_theta = Fc_phi * phi_theta
    
    Fq = -0.25 * 0.25 * o3 * (r * (4 * b * (2 * b2 + 9 * r2) * sphi &
       - 4 * r * (3 * (2 * b2 + r2) * phi &
       + br * Sin(3 * phi)) + r2 * r * Sin(4 * phi)))

    Fq_phi = - o3 * 0.25 * (r * (b * (2 * b2 + 9 * r2) * cphi &
        - 3 * r * (2 * b2 + r2 + br * Cos(3 * phi)) + r2 * r * Cos(4 * phi)))

    Fq_r = ( - (b * (2 * b2 + 27 * r2) * sphi) &
         + r * (12 * (b2 + r2) * phi + 3 * br * Sin(3 * phi) &
         - r2 * Sin(4 * phi))) * 0.25 * o3

    Fq_b = b * (r2 * phi - br * sphi * 0.5) &
        + r2 * 0.25 * r * (Sin(3 * phi) * o3 - 3 * sphi)
                
    Fq_rp = Fq_phi * phi_rp
    Fq_rm = Fq_phi * phi_rm
    Fq_bp = Fq_phi * phi_bp
    Fq_bpm = Fq_phi * phi_bpm
    Fq_theta = Fq_phi * phi_theta
        
    if (-c1 .eq. 2 * c2) then
        Fl = 0.d0
        Fl_phi = 0.d0
        Fl_r = 0.d0
        Fl_b = 0.d0        
    else          
        if (bpr .gt. 1.d0) then
        
            y = Sqrt(br)
            oy = 1.d0 / Sqrt(br)
            ox = 1.d0 / (b2 + r2 - 2 * br * cphi)
            
            alpha = 2 * y * (7 * r2 + b2 - 4.d0) * o9
            beta = -(3.d0 + 2*r*(b2 * b + 5 * b2 * r + 3*r*(-2.d0 + r2) + b*(-4.d0 + 7*r2))) * o9 * 0.5 * oy
            gamma = bpr * o3 * 0.5 * oy * obmr
            
            m = (1.d0 - bmr) * (1.d0 + bmr) * 0.25 * oy * oy
            n = ((bmr + 1.d0) * (bmr - 1.d0)) * obmr * obmr
            sqomm = Sqrt(1.d0 - m)

            ur = 4 * r * y
            vr = - r * (bpr + 1.d0) * (bpr - 1.d0) * oy
            ub = 2 * r * (b2 + r2 - 1.d0) * o3 * oy
            vb = vr * o3
            
            if (limbflag) then
            
                d = phi * o3 * 0.5 - Atan(bpr * tphihalf * obmr) * o3
                d_phi = (r2 - br * cphi) * o3 * ox
                d_r = - b * sphi * o3 * ox
                d_b = r * sphi * o3 * ox
                
                Fl_phi = 0.d0
                pr = 0.d0
                pb = 0.d0
                
                o = 1.d0
                ellippi = cel((sqomm), 1.d0 - n, (o), (o))
                ellipe = cel((sqomm), (o), (o), (1.d0 - m))
                ellipf = cel((sqomm), (o), (o), (o))
                
                eplusf = alpha * ellipe + beta * ellipf
                eplusf_r = ur * ellipe + vr * ellipf
                eplusf_b = ub * ellipe + vb * ellipf
            else                
                z = Sqrt((1.d0 - b) * (1.d0 + b) - r2 + 2 * br * cphi)
                oz = 1.d0 / z
                d = o3 * (phi * 0.5 - Atan(bpr * tphihalf * obmr)) &
                    - (2 * br * o9) * sphi * z
                d_phi = o3 * o3 * r * (3 * (r - b * cphi) * ox &
                    - 2 * b * cphi * z + 2 * b2 * r * sphi * sphi * oz)
                d_r = o9 * b * (-3.d0 * ox + 2 * (r2 - br * cphi) * oz - 2 * z) * sphi
                d_b = o9 * r * (3.d0 * ox + 2 * (b2 - br * cphi) * oz - 2 * z) * sphi
                
                Fl_phi = cphihalf * (-(n * (alpha + 4 * beta)) &
                       + 4 * m * (alpha + 2 * (beta + gamma)) &
                       + 4 * (m * alpha + n * beta) * cphi + n * alpha * Cos(2 * phi)) &
                       / (4 *  Sqrt((1.d0 + cphi) * (-1.d0 + 2 * m + cphi)) * (2 * m - n + n * cphi))
                    
                pr = b * sphi * (3.d0 - 4 * b2 + b2 * b2 - r2 * (4.d0 + r2) + 2 * br * (4.d0 - b2 + r2) * cphi) &
                   / (9 * y * Sqrt(2 * cphi - (b2 + r2 - 1.d0) * oy * oy) * (b2 + r2 - 2 * br * cphi))
                pb = r * sphi * (-3.d0 + 2 * b2 - b2 * b2 + 2 * r2 + r2 * r2 + 2 * br * (-2.d0 + b2 - r2) * cphi) &
                   / (9 * y * Sqrt(2 * cphi - (b2 + r2 - 1.d0) * oy * oy) * (b2 + r2 - 2 * br * cphi))
                                        
                tans = 1.d0 / Sqrt(m / (sphihalf * sphihalf) - 1.d0)
                
                o = 1.d0
                ellippi = el3((tans), (sqomm), 1.d0 - n)
                ellipe = el2((tans), (sqomm), (o), (1.d0 - m))
                ellipf = el2((tans), (sqomm), (o), (o))
                
                eplusf = alpha * ellipe + beta * ellipf
                eplusf_r = ur * ellipe + vr * ellipf
                eplusf_b = ub * ellipe + vb * ellipf

            end if
            Fl = eplusf + gamma * ellippi + d
            Fl_phi = Fl_phi + d_phi  
            Fl_r = eplusf_r + pr + d_r
            Fl_b = eplusf_b + pb + d_b
            
        else
        
            z = Sqrt((1.d0 - bmr) * (1.d0 + bmr))
            oz = 1.d0 / z
            y = Sqrt((1.d0 - b) * (1.d0 + b) - r2 + 2 * br * cphi)
            oy = 1.d0 / y
            x = (b2 + r2 - 2 * br * cphi)
            ox = 1.d0 / x
            
            alpha = (7 * r2 + b2 - 4.d0) * z * o9
            beta = (r2 * r2 + b2 * b2 + r2 - b2 * (5.d0 + 2 * r2) + 1.d0) * o9 * oz
            gamma = bpr * o3 * oz * obmr
                        
            n = -4 * br * obmr * obmr
            m = 4 * br * oz * oz
            
            sqomm = Sqrt(1.d0 - m)
            
            ur = 2 * r * z
            pr = b * sphi * (3.d0 - 4 * b2 + b2 * b2 - r2 * (4.d0 + r2) + 2 * br * (4.d0 - b2 + r2) * cphi) &
               / (9 * y * (1.d0 - y) * (1.d0 + y))
                
            ub = z * ((b + 1.d0) * (b - 1.d0) + r2) / (3 * b)
            vb = (b2 * b2 + ((r - 1.d0) * (r + 1.d0))**2.d0 - 2 * b2 * (1.d0 + r2)) / (3 * b * z)
            pb = r * sphi * (-3.d0 + 2 * b2 - b2 * b2 + 2 * r2 + r2 * r2 &
               + 2 * br * (-2.d0 + b2 - r2) * cphi) * ox * o9 * oy 
            
            d = o3 * (phi * 0.5 - Atan(bpr * tphihalf * obmr)) &
                - 2 * br * o9 * sphi * y
            d_r = o9 * b * sphi * (-3.d0 * ox + 2 * (r2 - br * cphi) * oy - 2 * y)
            d_b = o9 * r * (3.d0 * ox + 2 * (b2 - br * cphi) * oy - 2 * y) * sphi
            d_phi = o9 * r * (3 * (r - b * cphi) * ox &
                - 2 * b * cphi * y + 2 * b2 * r * sphi * sphi * oy)
             
            o = 1.d0
            ellipe = el2((tphihalf), (sqomm), (o), (1.d0 - m))
            ellipf = el2((tphihalf), (sqomm), (o), (o))
            
            if (b .eq. r) then
                Fl_phi = ((2.d0 + m * (cphi - 1.d0)) * alpha + 2 * beta) &
                       / (2 * Sqrt(4.d0 + 2 * m * (cphi - 1.d0)))
                ellippi = 0.d0
                gamma = 0.d0
            else
                Fl_phi = ((2 - n + m * (cphi - 1.d0)) * alpha + (2 - n) * beta + 2 * gamma &
                    + n * cphi * (alpha + beta) + 2 * m * n * alpha * sphihalf**4) &
                    / (Sqrt(2 * (2.d0 - m + m * cphi)) * (2.d0 - n + n * cphi))
                ellippi = el3((tphihalf), (sqomm), (1.d0 - n))
            end if
                
            eplusf = alpha * ellipe + beta * ellipf
            eplusf_r = ur * ellipe
            eplusf_b = ub * ellipe + vb * ellipf 
            
            Fl = eplusf + gamma * ellippi + d
            Fl_phi = Fl_phi + d_phi
            Fl_r = eplusf_r + pr + d_r
            Fl_b = eplusf_b + pb + d_b
        end if
        
        if (b .eq. 0.d0) then
            Fl = phi * (1.d0 - ((1.d0 + r) * (1.d0 - r)) ** (1.5)) * o3
            Fl_phi = (1.d0 - ((1.d0 + r) * (1.d0 - r)) ** (1.5)) * o3
            Fl_r = -pi * r * Sqrt((1.d0 + r) * (1.d0 - r))
        end if
        
    end if
    
2   cc = 1.d0 - c1 - 2 * c2
    cl = c1 + 2 * c2
    cq = c2
        
    Fl_rp = Fl_phi * phi_rp
    Fl_rm = Fl_phi * phi_rm
    Fl_bp = Fl_phi * phi_bp
    Fl_bpm = Fl_phi * phi_bpm
    Fl_theta = Fl_phi * phi_theta
    
    if (pflag) then
        Fq_bp = Fq_bp + Fq_b
        Fq_rp = Fq_rp + Fq_r
        
        Fc_bp = Fc_bp + Fc_b
        Fc_rp = Fc_rp + Fc_r
        
        Fl_bp = Fl_bp + Fl_b
        Fl_rp = Fl_rp + Fl_r
    else
        Fq_theta = Fq_theta + Fq_b * bm_theta
        Fq_bpm = Fq_bpm + Fq_b * bm_bpm
        Fq_bp = Fq_bp + Fq_b * bm_bp
        Fq_rm = Fq_rm + Fq_r
        
        Fc_theta = Fc_theta + Fc_b * bm_theta
        Fc_bpm = Fc_bpm + Fc_b * bm_bpm
        Fc_bp = Fc_bp + Fc_b * bm_bp
        Fc_rm = Fc_rm + Fc_r
        
        Fl_theta = Fl_theta + Fl_b * bm_theta
        Fl_bpm = Fl_bpm + Fl_b * bm_bpm
        Fl_bp = Fl_bp + Fl_b * bm_bp
        Fl_rm = Fl_rm + Fl_r
    end if
    
    F(1) = cc * Fc + cl * Fl + cq * Fq
    F(2) = cc * Fc_rp + cl * Fl_rp + cq * Fq_rp
    F(3) = cc * Fc_rm + cl * Fl_rm + cq * Fq_rm
    F(4) = cc * Fc_bp + cl * Fl_bp + cq * Fq_bp
    F(5) = cc * Fc_bpm + cl * Fl_bpm + cq * Fq_bpm
    F(6) = cc * Fc_theta + cl * Fl_theta + cq * Fq_theta
    F(7) = - Fc + Fl 
    F(8) = -2 * Fc + 2 * Fl + Fq

    return
end function

end module phot