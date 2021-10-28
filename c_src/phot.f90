module phot
use iso_c_binding
use ellip

implicit none

real*8, parameter :: pi = 3.14159265358979323846, pihalf = 1.5707963267948966, twopithree = 2.0943951023931953
real*8, parameter :: o3 = 0.33333333333333333333, o9 = 0.1111111111111111111, twopi = 6.283185307179586
real*8, parameter :: pithird = 1.0471975511965976, pisixth = 0.5235987755982988

contains

! computes the angular location of the intersections between moon and planet from 
! from the perspective of each body 
subroutine compute_phis(rp, rm, bp, bm, bpm, phim1, phim2, phip1, phip2, &
    phim1_bpm, phim1_bp, phim1_bm, phim1_rp, phim1_rm, phim2_bpm, phim2_bp, &
    phim2_bm, phim2_rp, phim2_rm, phip1_bpm, phip1_bp, phip1_bm, &
    phip1_rp, phip1_rm, phip2_bpm, phip2_bp, phip2_bm, phip2_rp, phip2_rm)

    real*8 :: rp, rm, rp2, rm2, bp, bm, bpm, bpm2, denom
    real*8, intent(out) :: phim1, phim2, phip1, phip2
    real*8, intent(out) :: phim1_bpm, phim1_rp, phim1_rm, phip1_bpm, phip1_rm, phip1_rp
    real*8, intent(out) :: phim2_bpm, phim2_rp, phim2_rm, phip2_bpm, phip2_rm, phip2_rp
    real*8, intent(out) :: phim1_bm, phim1_bp, phim2_bm, phim2_bp, phip1_bm, phip1_bp, phip2_bm, phip2_bp
    real*8 :: a, b, c, thetam, thetap, tmp, area, phip, phim
    real*8 :: thetam_bpm, thetap_bpm, thetam_rm, thetap_rm, thetam_rp, thetap_rp
    real*8 :: phim_bpm, phip_bpm, x, y, phim_bp, phim_bm, phip_bp, phip_bm

    bpm2 = bpm * bpm
    rp2 = rp * rp
    rm2 = rm * rm
    
    ! first compute the angle from the planet center - moon center line 
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
    
    ! then compute the angle to the planet or moon from the x-axis 
    ! in order to transform into the correct coordinates 
    
    ! need a fix for when planet/moon are perfectly aligned. 
    ! The angles are right but the derivatives aren't -- they blow up 
    ! at bpm = bp - bm and bpm = bm - bp
    if (bpm .eq. bm - bp) then
        phim = 0.d0
        phim_bpm = 0.d0
        phim_bm = 0.d0
        phim_bp = 0.d0
        
        phip = pi
        phip_bpm = 0.d0
        phip_bm = 0.d0
        phim_bp = 0.d0
    else if (bpm .eq. bp - bm) then
        phim = 0.d0
        phim_bpm = 0.d0
        phim_bm = 0.d0
        phim_bp = 0.d0
        
        phip = pi
        phip_bpm = 0.d0
        phip_bm = 0.d0
        phim_bp = 0.d0
    else
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
        phim_bm = ((bpm - bm) * (bpm + bm) - bp * bp) / (bm * area)
        phim_bp = 2 * bp / area
    
        phip = Atan2(area, (bp - bm) * (bp + bm) + bpm2)
        phip_bpm = ((bp - bm) * (bp + bm) - bpm2) / (bpm * area)
        phip_bm = 2 * bm / area
        phip_bp = ((bpm + bp) * (bpm - bp) - bm * bm) / (bp * area)
    end if

    phim1 = phim + thetam
    phim1_bpm = phim_bpm + thetam_bpm
    phim1_bp = phim_bp
    phim1_bm = phim_bm
    phim1_rp = thetam_rp
    phim1_rm = thetam_rm
    
    phim2 = phim - thetam
    phim2_bpm = phim_bpm - thetam_bpm
    phim2_bm = phim_bm
    phim2_bp = phim_bp
    phim2_rp = -thetam_rp
    phim2_rm = -thetam_rm
    
    phip1 = phip + thetap
    phip1_bpm = phip_bpm + thetap_bpm
    phip1_bp = phip_bp
    phip1_bm = phip_bm
    phip1_rp = thetap_rp
    phip1_rm = thetap_rm
    
    phip2 = phip - thetap
    phip2_bpm = phip_bpm - thetap_bpm
    phip2_bp = phip_bp
    phip2_bm = phip_bm
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

! computes the angle to the moon-star or planet-star intersections from the 
! perspective of the moon/planet and the star
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

! main loop to compute the flux at each timestep by finding the correct geometry and
! calling the integration routines 
subroutine flux(c1, c2, rp, rm, bp2, bm2, bpm2, lc, j) bind(C, name="flux")

    integer (c_int), bind(C) :: j
    integer :: i
    real (c_double), bind(C) :: rp, rm
    real (c_double), bind(C), dimension(j) :: bp2, bm2, bpm2
    real (c_double), bind(C), intent(out), dimension(8, j) :: lc
    real*8, dimension(8) :: f0
    real (c_double), bind(C) :: c1, c2
    
    real*8 :: theta, thetapm, phim1, phim2, phip1, phip2, thetap, thetam
    real*8 :: phim1_bpm, phim1_bp, phim1_bm, phim1_rp, phim1_rm, phip1_bpm, phip1_bp, phip1_bm, phip1_rm, phip1_rp
    real*8 :: phim2_bpm, phim2_bp, phim2_bm, phip2_bp, phip2_bm, phim2_rp, phim2_rm, phip2_bpm, phip2_rm, phip2_rp
    real*8 :: phim, phip, cosphim, cosphip, costhetapm, costheta, cosphi1, cosphi2
    real*8 :: bp2i, bm2i, bpm2i, bpi, bmi, bpmi, rp2, rm2, obpi, obmi
    real*8 :: thetap_bp, thetap_rp, phip_bp, phip_rp, thetam_bp, thetam_rp, phim_bp, phim_rp
    real*8 :: theta_bp, theta_rp, phi_bp, phi_rp, phi, thetam_bm, thetam_rm, phi_bm, phi_rm
    real*8 :: d1, d2, area, phim_bm, phim_rm, theta_bm, theta_rm, theta_bpm, thetapm_bpm
    real*8 :: of0, a, b, c, tmp, thetapm_bp, thetapm_bm
    
    real*8 :: bp(j), bm(j), bpm(j)
    bp = Sqrt(bp2)
    bm = Sqrt(bm2)
    bpm = Sqrt(bpm2)
    
    rp2 = rp * rp
    rm2 = rm * rm
    
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
    
        bp2i = bp2(i)
        bm2i = bm2(i)
        bpm2i = bpm2(i)
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
                    if (bmi + rm .le. 1.d0) then
                        ! moon completely overlaps star, planet is outside of star
                        lc(:, i) = (f0 - 2 * Fcomplete(c1, c2, rm, bmi, .FALSE.)) * of0
                    else
                        ! moon partially overlaps star, planet is outside of star
                        call compute_theta(rm, bmi, theta, phi, theta_bm, theta_rm, phi_bm, phi_rm)
                        lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, 0.d0, 0.d0, -phi_bm, -phi_rm, 0.d0) &
                                - F(c1, c2, theta, rm, bmi, 0.d0, 0.d0, theta_bm, theta_rm, 0.d0, .FALSE.)) * of0
                    end if
                end if
            else
                if (bmi .gt. rm + 1.d0) then
                    if (bpi + rp .le. 1.d0) then
                        ! planet completely overlaps star, moon is outside of star
                        lc(:, i) = (f0 - 2 * Fcomplete(c1, c2, rp, bpi, .TRUE.)) * of0
                    else
                        ! planet partially overlaps star, moon is outside of star
                        call compute_theta(rp, bpi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                        lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, -phi_bp, -phi_rp, 0.d0, 0.d0, 0.d0) &
                              - F(c1, c2, theta, rp, bpi, theta_bp, theta_rp, 0.d0, 0.d0, 0.d0, .TRUE.)) * of0
                    end if
                else
                    if (bpi + rp .le. 1.d0) then
                        if (bmi + rm .le. 1.d0) then
                            ! moon and planet both completely overlap star, they do not overlap each othe
                            lc(:, i) = (f0 - 2 * (Fcomplete(c1, c2, rm, bmi, .FALSE.) &
                                  + Fcomplete(c1, c2, rp, bpi, .TRUE.))) * of0
                        else
                            ! planet completely overlaps star, moon partially overlaps star, they do not overlap each other
                            call compute_theta(rm, bmi, theta, phi, theta_bm, theta_rm, phi_bm, phi_rm)
                            lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, 0.d0, 0.d0, -phi_bm, -phi_rm, 0.d0) &
                                  - F(c1, c2, theta, rm, bmi, 0.d0, 0.d0, theta_bm, theta_rm, 0.d0, .FALSE.) &
                                  - Fcomplete(c1, c2, rp, bpi, .TRUE.)) * of0
                        end if
                    else
                        if (bmi + rm .le. 1.d0) then
                            ! planet partially overlaps star, moon fully overlaps star, they do not overlap each other
                            call compute_theta(rp, bpi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                            lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, -phi_bp, -phi_rp, 0.d0, 0.d0, 0.d0) &
                                  - F(c1, c2, theta, rp, bpi, theta_bp, theta_rp, 0.d0, 0.d0, 0.d0, .TRUE.) & 
                                  - Fcomplete(c1, c2, rm, bmi, .FALSE.)) * of0
                        else
                            ! moon and planet both partially overlap star, but not each other
                            call compute_theta(rp, bpi, thetap, phip, thetap_bp, thetap_rp, phip_bp, phip_rp)
                            call compute_theta(rm, bmi, thetam, phim, thetam_bm, thetam_rm, phim_bm, phim_rm)
                            phi = phip + phim
                            lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, -phip_bp, -phip_rp, -phim_bm, -phim_rm, 0.d0) &
                                  - F(c1, c2, thetap, rp, bpi, thetap_bp, thetap_rp, 0.d0, 0.d0, 0.d0, .TRUE.) & 
                                  - F(c1, c2, thetam, rm, bmi, 0.d0, 0.d0, thetam_bm, thetam_rm, 0.d0, .FALSE.)) * of0
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
                    call compute_theta(rm, bmi, theta, phi, theta_bm, theta_rm, phi_bm, phi_rm)
                    lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, 0.d0, 0.d0, -phi_bm, -phi_rm, 0.d0) &
                            - F(c1, c2, theta, rm, bmi, 0.d0, 0.d0, theta_bm, theta_rm, 0.d0, .FALSE.)) * of0
                end if
            else
                if (bmi .gt. rm + 1.d0) then
                    if (bpi + rp .le. 1.d0) then
                        ! planet fully overlaps star, moon does not overlap star
                        lc(:, i) = (f0 - 2 * Fcomplete(c1, c2, rp, bpi, .TRUE.)) * of0
                    else
                        ! planet partially overlaps star, moon does not overlap star
                        call compute_theta(rp,  bpi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                        lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, -phi_bp, -phi_rp, 0.d0, 0.d0, 0.d0) &
                              - F(c1, c2, theta, rp, bpi, theta_bp, theta_rp, 0.d0, 0.d0, 0.d0, .TRUE.)) * of0
                    end if
                else
                    if (bpi + rp .le. 1.d0) then
                        if (bmi + rm .le. 1.d0) then
                            if (bpmi + rm .le. rp) then
                                ! moon and planet both overlap star, moon fully overlapped by planet
                                lc(:, i) = (f0 - 2 * Fcomplete(c1, c2, rp, bpi, .TRUE.)) * of0
                            else
                                ! moon and planet both overlap star, moon and planet partially overlap each other 
                                call compute_phis(rp, rm, bpi, bmi, bpmi, phim1, phim2, phip1, phip2, &
                                    phim1_bpm, phim1_bp, phim1_bm, phim1_rp, phim1_rm, &
                                    phim2_bpm, phim2_bp, phim2_bm, phim2_rp, phim2_rm, &
                                    phip1_bpm, phip1_bp, phip1_bm, phip1_rp, phip1_rm, &
                                    phip2_bpm, phip2_bp, phip2_bm, phip2_rp, phip2_rm)
                                lc(:, i) = (f0 - (Arc(c1, c2, phip1, phip2, rp, bpi, &
                                                    phip1_bp, phip1_rp, phip1_bm, phip1_rm, phip1_bpm, &
                                                    phip2_bp, phip2_rp, phip2_bm, phip2_rm, phip2_bpm, .TRUE.) &
                                                + Arc(c1, c2, phim1, phim2, rm, bmi, &
                                                      phim1_bp, phim1_rp, phim1_bm, phim1_rm, phim1_bpm, &
                                                      phim2_bp, phim2_rp, phim2_bm, phim2_rm, phim2_bpm, .False.))) * of0
                            end if
                        else 
                            if (bpmi + rm .le. rp) then
                                ! not physical (planet fully overlaps star, moon partially 
                                ! overlaps star while fully overlapped by planet.)
                            else
                                ! planet fully overlaps star, moon partially overlaps star and partially overlaps planet. 
                                call compute_phis(rp, rm, bpi, bmi, bpmi, phim1, phim2, phip1, phip2, &
                                    phim1_bpm, phim1_bp, phim1_bm, phim1_rp, phim1_rm, &
                                    phim2_bpm, phim2_bp, phim2_bm, phim2_rp, phim2_rm, &
                                    phip1_bpm, phip1_bp, phip1_bm, phip1_rp, phip1_rm, &
                                    phip2_bpm, phip2_bp, phip2_bm, phip2_rp, phip2_rm)
                                call compute_theta(rm, bmi, theta, phi, theta_bm, theta_rm, phi_bm, phi_rm)
                                lc(:, i) = (2 * Fstar(c1, c2, pi - phi, 0.d0, 0.d0, -phi_bm, -phi_rm, 0.d0) &
                                  - Arc(c1, c2, -theta, phim2, rm, bmi, &
                                        0.d0, 0.d0, -theta_bm, -theta_rm, 0.d0, &
                                        phim2_bp, phim2_rp, phim2_bm, phim2_rm, phim2_bpm, .FALSE.) &
                                  - Arc(c1, c2, phim1, theta, rm, bmi, &
                                        phim2_bp, phim1_rp, phim1_bm, phim1_rm, phim1_bpm, &
                                        0.d0, 0.d0, theta_bm, theta_rm, 0.d0, .FALSE.) &
                                  - Arc(c1, c2, phip1, phip2, rp, bpi, &
                                        phip1_bp, phip1_rp, phip1_bm, phip1_rm, phip1_bpm, &
                                        phip2_bp, phip2_rp, phip2_bm, phip2_rm, phip2_bpm, .TRUE.)) * of0
                            end if
                        end if
                    else
                        if (bmi + rm .le. 1.d0) then 
                            if (bpmi + rm .le. rp) then
                                ! planet partially overlaps star, moon fully overlaps star but is completely overlapped by planet 
                                call compute_theta(rp, bpi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                                lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, -phi_bp, -phi_rp, 0.d0, 0.d0, 0.d0) &
                                    - F(c1, c2, theta, rp, bpi, theta_bp, theta_rp, 0.d0, 0.d0, 0.d0, .TRUE.)) * of0
                            else
                                ! planet partially overlaps star, moon fully overlaps star and only partially overlaps planet
                                call compute_phis(rp, rm, bpi, bmi, bpmi, phim1, phim2, phip1, phip2, &
                                    phim1_bpm, phim1_bp, phim1_bm, phim1_rp, phim1_rm, &
                                    phim2_bpm, phim2_bp, phim2_bm, phim2_rp, phim2_rm, &
                                    phip1_bpm, phip1_bp, phip1_bm, phip1_rp, phip1_rm, &
                                    phip2_bpm, phip2_bp, phip2_bm, phip2_rp, phip2_rm)
                                call compute_theta(rp, bpi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                                lc(:, i) = (2 * Fstar(c1, c2, pi - phi, -phi_bp, -phi_rp, 0.d0, 0.d0, 0.d0) &
                                  - Arc(c1, c2, -theta, phip2, rp, bpi, &
                                        -theta_bp, -theta_rp, 0.d0, 0.d0, 0.d0, &
                                        phip2_bp, phip2_rp, phip2_bm, phip2_rm, phip2_bpm, .TRUE.) &
                                  - Arc(c1, c2, phip1, theta, rp, bpi, &
                                        phip1_bp, phip1_rp, phip1_bm, phip1_rm, phip1_bpm, &
                                        theta_bp, theta_rp, 0.d0, 0.d0, 0.d0, .TRUE.) &
                                  - Arc(c1, c2, phim1, phim2, rm, bmi, & 
                                        phim1_bp, phim1_rp, phim1_bm, phim1_rm, phim1_bpm, &
                                        phim2_bp, phim2_rp, phim2_bm, phim2_rm, phim2_bpm, .FALSE.)) * of0
                            end if
                        else
                            obpi = 1.d0 / bpi
                            obmi = 1.d0 / bmi
                            if (bpmi + rm .le. rp) then
                                ! planet and moon both partially overlap star but moon is fully overlapped by the planet
                                call compute_theta(rp, bpi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                                lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, -phi_bp, -phi_rp, 0.d0, 0.d0, 0.d0) &
                                    - F(c1, c2, theta, rp, bpi, theta_bp, theta_rp, 0.d0, 0.d0, 0.d0, .TRUE.)) * of0
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
                                thetapm_bm = ((bpi + bmi) * (bpi - bmi) - bpmi * bpmi) / (bmi * area)
                                thetapm_bp = ((bmi + bpi) * (bmi - bpi) - bpmi * bpmi) / (bpi * area)
                                
                                if (thetapm + phim .le. phip) then
                                        call compute_phis(rp, rm, bpi, bmi, bpmi, phim1, phim2, phip1, phip2, &
                                                phim1_bpm, phim1_bp, phim1_bm, phim1_rp, phim1_rm, &
                                                phim2_bpm, phim2_bp, phim2_bm, phim2_rp, phim2_rm, &
                                                phip1_bpm, phip1_bp, phip1_bm, phip1_rp, phip1_rm, &
                                                phip2_bpm, phip2_bp, phip2_bm, phip2_rp, phip2_rm)
                                        call compute_theta(rp,  bpi, theta, phi, theta_bp, theta_rp, phi_bp, phi_rp)
                                        if (phip2 .gt. theta) then
                                            ! planet and moon both partially overlap the star and each other but the 
                                            ! moon-star overlap is contained within the planet-star overlap
                                            lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, -phi_bp, -phi_rp, 0.d0, 0.d0, 0.d0) &
                                                - F(c1, c2, theta, rp, bpi, theta_bp, theta_rp, 0.d0, 0.d0, 0.d0, .TRUE.)) * of0
                                        else
                                            ! planet and moon both partially overlap star and each other but the 
                                            ! planet-star intersections are overlapped by the planet
                                            lc(:, i) = (2 * Fstar(c1, c2, pi - phi, -phi_bp, -phi_rp, 0.d0, 0.d0, 0.d0) &
                                              - Arc(c1, c2, -theta, phip2, rp, bpi, &
                                                    -theta_bp, -theta_rp, 0.d0, 0.d0, 0.d0, &
                                                    phip2_bp, phip2_rp, phip2_bm, phip2_rm, phip2_bpm, .TRUE.) &
                                              - Arc(c1, c2, phip1, theta, rp, bpi, &
                                                    phip1_bp, phip1_rp, phip1_bm, phip1_rm, phip1_bpm, &
                                                    theta_bp, theta_rp, 0.d0, 0.d0, 0.d0, .TRUE.) &
                                              - Arc(c1, c2, phim1, phim2, rm, bmi, &
                                                    phim1_bp, phim1_rp, phim1_bm, phim1_rm, phim1_bpm, &
                                                    phim2_bp, phim2_rp, phim2_bm, phim2_rm, phim2_bpm, .FALSE.)) * of0
                                        end if
                                else if (thetapm + phip .le. phim) then
                                    if ((bpi - rp) .le. (bmi - rm)) then
                                        ! planet and moon both partially overlap the star and each other but the 
                                        ! planet-star intersections are overlapped by the moon
                                        ! I'm not sure this is physical either -- can you draw a diagram where 
                                        ! the moon overlaps both of the planet-star 
                                        ! intersections without the planet-star overlap being entirely within 
                                        ! the moon-star region of overlap?
                                        call compute_phis(rp, rm, bpi, bmi, bpmi, phim1, phim2, phip1, phip2, &
                                                phim1_bpm, phim1_bp, phim1_bm, phim1_rp, phim1_rm, &
                                                phim2_bpm, phim2_bp, phim2_bm, phim2_rp, phim2_rm, &
                                                phip1_bpm, phip1_bp, phip1_bm, phip1_rp, phip1_rm, &
                                                phip2_bpm, phip2_bp, phip2_bm, phip2_rp, phip2_rm)
                                        call compute_theta(rm, bmi, theta, phi, theta_bm, theta_rm, phi_bm, phi_rm)
                                        lc(:, i) = (2 * Fstar(c1, c2, pi - phi, 0.d0, 0.d0, -phi_bm, -phi_rm, 0.d0) &
                                            - Arc(c1, c2, -theta, phim2, rm, bmi, &
                                                  0.d0, 0.d0, -theta_bm, -theta_rm, 0.d0, &
                                                  phim2_bp, phim2_rp, phim2_bm, phim2_rm, phim2_bpm, .FALSE.) &
                                            - Arc(c1, c2, phim1, theta, rm, bmi, &
                                                  phim1_bp, phim1_rp, phim1_bm, phim1_rm, phim1_bpm, &
                                                  0.d0, 0.d0, theta_bm, theta_rm, 0.d0, .FALSE.) & 
                                            - Arc(c1, c2, phip1, phip2, rp, bpi, &
                                                  phip1_bp, phip1_rp, phip1_bm, phip1_rm, phip1_bpm, &
                                                  phip2_bp, phip2_rp, phip2_bm, phip2_rm, phip2_bpm, .TRUE.)) * of0
                                    else
                                        ! planet and moon both partially overlap the star and each other but 
                                        ! the planet-star overlap is  entirely within the moon-star overlap
                                        call compute_theta(rm, bmi, theta, phi, theta_bm, theta_rm, phi_bm, phi_rm)
                                        lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, 0.d0, 0.d0, -phi_bm, -phi_rm, 0.d0) &
                                          - F(c1, c2, theta, rm, bmi, 0.d0, 0.d0, theta_bm, theta_rm, 0.d0, .FALSE.)) * of0
                                    end if
                                else
                                    costheta = (bpm2i + bm2i - bp2i) / (2 * bpmi * bmi)
                                    cosphim = (bpm2i + rm2 - rp2) / (2 * bpmi * rm)
                                    cosphi1 = Cos(Acos(costheta) - Acos(cosphim))
                                    cosphi2 = Cos(Acos(costheta) + Acos(cosphim))
                                    d1 = rm2 + bm2i - 2 * rm * bmi * cosphi1
                                    d2 = rm2 + bm2i - 2 * rm * bmi * cosphi2
                                    if (d1 .gt. 1.d0) then
                                        ! planet and moon both partially overlap star and each other, 
                                        ! but the planet/moon overlap does not overlap the star
                                        call compute_theta(rp, bpi, thetap, phip, thetap_bp, thetap_rp, phip_bp, phip_rp)
                                        call compute_theta(rm, bmi, thetam, phim, thetam_bm, thetam_rm, phim_bm, phim_rm)
                                        phi = phip + phim
                                        lc(:, i) = 2 * (Fstar(c1, c2, pi - phi, -phip_bp, -phip_rp, -phim_bm, -phim_rm, 0.d0) &
                                          - F(c1, c2, thetap, rp, bpi, theta_bp, thetap_rp, 0.d0, 0.d0, 0.d0, .TRUE.) & 
                                          - F(c1, c2, thetam, rm, bmi, 0.d0, 0.d0, thetam_rm, thetam_rp, 0.d0, .FALSE.)) * of0
                                    else if (d2 .le. 1.d0) then
                                        ! planet and moon both partially overlap star and each other, 
                                        ! with the planet/moon overlap fully overlapping the star
                                        call compute_phis(rp, rm, bpi, bmi, bpmi, phim1, phim2, phip1, phip2, &
                                                phim1_bpm, phim1_bp, phim1_bm, phim1_rp, phim1_rm, &
                                                phim2_bpm, phim2_bp, phim2_bm, phim2_rp, phim2_rm, &
                                                phip1_bpm, phip1_bp, phip1_bm, phip1_rp, phip1_rm, &
                                                phip2_bpm, phip2_bp, phip2_bm, phip2_rp, phip2_rm)
                                        call compute_theta(rm, bmi, thetam, phim, thetam_bm, thetam_rm, phim_bm, phim_rm)
                                        call compute_theta(rp, bpi, thetap, phip, thetap_bp, thetap_rp, phip_bp, phip_rp)
                                        phi = phip + phim
                                        lc(:, i) = (2 * Fstar(c1, c2, pi - phi, -phip_bp, -phip_rp, -phim_bm, -phim_rm, 0.d0) &
                                            - Arc(c1, c2, -thetam, -phim1, rm, bmi, &
                                                  0.d0, 0.d0, -thetam_bm, -thetam_rm, 0.d0, &
                                                  -phim1_bp, -phim1_rp, -phim1_bm, -phim1_rm, -phim1_bpm, .FALSE.) & 
                                            - Arc(c1, c2, -phim2, thetam, rm, bmi, &
                                                  -phim2_bp, -phim2_rp, -phim2_bm, -phim2_rm, -phim2_bpm, &
                                                  0.d0, 0.d0, theta_bm, theta_rm, 0.d0, .FALSE.) &
                                            - Arc(c1, c2, phip1, thetap, rp, bpi, &
                                                  phip1_bp, phip1_rp, phip1_bm, phip1_rm, phip1_bpm, &
                                                  thetap_bp, thetap_rp, 0.d0, 0.d0, 0.d0, .TRUE.) &
                                            - Arc(c1, c2, -thetap, phip2, rp, bpi, &
                                                  -thetap_bp, -thetap_rp, 0.d0, 0.d0, 0.d0, &
                                                  phip2_bp, phip2_rp, phip2_bm, phip2_rm, phip2_bpm, .TRUE.)) * of0
                                    else
                                        ! planet and moon both partially overlap star and each other, 
                                        ! with the planet/moon overlap partially overlapping the star
                                        call compute_phis(rp, rm, bpi, bmi, bpmi, phim1, phim2, phip1, phip2, &
                                                phim1_bpm, phim1_bp, phim1_bm, phim1_rp, phim1_rm, &
                                                phim2_bpm, phim2_bp, phim2_bm, phim2_rp, phim2_rm, &
                                                phip1_bpm, phip1_bp, phip1_bm, phip1_rp, phip1_rm, &
                                                phip2_bpm, phip2_bp, phip2_bm, phip2_rp, phip2_rm)
                                        call compute_theta(rm, bmi, thetam, phim, thetam_bm, thetam_rm, phim_bm, phim_rm)
                                        call compute_theta(rp, bpi, thetap, phip, thetap_bp, thetap_rp, phip_bp, phip_rp)
                                        phi = 0.5 * (phip + phim + thetapm)
                                        lc(:, i) = (2 * Fstar(c1, c2, pi - phi, -0.5 * (phip_bp + thetapm_bp), &
                                                           -0.5 * phip_rp, -0.5 * (phim_bm + thetapm_bm), &
                                                           -0.5 * phim_rm, -0.5 * thetapm_bpm) & 
                                            - Arc(c1, c2, -phim2, thetam, rm, bmi, &
                                                  -phim2_bp, -phim2_rp, -phim2_bm, -phim2_rm, -phim2_bpm, &
                                                  0.d0, 0.d0, thetam_bm, thetam_rm, 0.d0, .FALSE.) & 
                                            - Arc(c1, c2, -thetap, phip2, rp, bpi, &
                                                  -thetap_bp, -thetap_rp, 0.d0, 0.d0, 0.d0, &
                                                  phip2_bp, phip2_rp, phip2_bm, phip2_rm, phip2_bpm, .TRUE.)) * of0
                                    end if
                                end if
                            end if
                        end if
                    end if
                end if
            end if  
        end if
        lc(:, i) = lc(:, i) - f0 * of0
    end do
    return
    
end

! work out the right sign and order of the integration and call the integration routine 
! to integrate along an arbitrary arc of the planet or moon 
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

! integrate along the limb of the star
function Fstar(c1, c2, phi, phi_bp, phi_rp, phi_bm, phi_rm, phi_bpm)

    real*8, dimension(8) :: Fstar

    real*8 :: c1, c2, phi, Fc, Fq, Fl
    real*8 :: Fc_bp, Fc_rp, Fc_bm, Fc_rm, Fc_bpm
    real*8 :: Fq_bp, Fq_rp, Fq_bm, Fq_rm, Fq_bpm
    real*8 :: Fl_bp, Fl_rp, Fl_bm, Fl_rm, Fl_bpm
    real*8 :: Fc_phi, Fq_phi, Fl_phi
    real*8 :: phi_bp, phi_rp, phi_bm, phi_rm, phi_bpm
    real*8 :: cc, cl, cq
    
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
    
    cc = 1.d0 - c1 - 2 * c2
    cl = c1 + 2 * c2
    cq = c2
        
    Fstar(1) = cc * Fc + cl * Fl + cq * Fq
    Fstar(2) = cc * Fc_bp + cl * Fl_bp + cq * Fq_bp
    Fstar(3) = cc * Fc_rp + cl * Fl_rp + cq * Fq_rp
    Fstar(4) = cc * Fc_bm + cl * Fl_bm + cq * Fq_bm
    Fstar(5) = cc * Fc_rm + cl * Fl_rm + cq * Fq_rm
    Fstar(6) = cc * Fc_bpm + cl * Fl_bpm + cq * Fq_bpm
    Fstar(7) = - Fc + Fl 
    Fstar(8) = -2 * Fc + 2 * Fl + Fq
    return
    
end function

! integrate around the entire planet/moon 
function Fcomplete(c1, c2, r, b, pflag)

    logical :: pflag
    real*8, dimension(8) :: Fcomplete
    
    real*8 :: c1, c2, r, b, Fc, Fq, Fl
    real*8 :: Fc_r, Fq_r, Fl_r, Fc_b, Fq_b, Fl_b
    real*8 :: o, sqomm, nmo, mmn, mmo, obmr
    real*8 :: cc, cl, cq
    real*8 :: gamma, n, m, x, y, tans, sphi, br, bmr, bpr
    real*8 :: r2, b2, apb, apbo
    real*8 :: ellippi, eplusf
    real*8 :: ellippi_r, eplusf_r, alpha, alpha_r, beta, beta_r
    real*8 :: ur, vr, pr, qr, m_r, n_r, gamma_r, denom
    real*8 :: ub, vb, pb, qb, m_b, n_b, gamma_b, alpha_b, beta_b, eplusf_b, ellippi_b
    real*8 :: ellipe, ellipf
    
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
        if (b == r) then
            if (r == 0.5) then
            
                Fl = pisixth + 2.d0 * o9
                Fl_r = 0.d0
                Fl_b = 0.d0                
            else
            
                ! much to fix here! shouldn't be any n's!!! 
                m = 4 * r2
                m_r = 8 * r
                sqomm = Sqrt(1.d0 - m)
                
                alpha = 4 * (2 * r2 - 1.d0) * o9
                alpha_r = 16 * r * o9
                
                beta = (1.d0 - 4 * r2) * o9
                beta_r = -8 * r * o9
                
                ur = ((m - n) * (1.d0 - n) * ((m - 1.d0) * (m_r * alpha + 2 * m * alpha_r) + m_r * beta) &
                    - m * (m_r - m_r * n)) &
                    / (2 * (m - 1.d0) * m * (m - n) * (1.d0 - n))
                vr = (2 * beta_r * m * (1.d0 - n) * n - m_r * (1.d0 - n) * n * (alpha + beta)) &
                    / (2 * m * (1.d0 - n) * n)
                
                eplusf = cel((sqomm), o, alpha + beta, alpha * (1.d0 - m) + beta)
                eplusf_r = cel((sqomm), (o), ur + vr, ur * (1.d0 - m) + vr)
                
                Fl = eplusf + pisixth
                Fl_r = eplusf_r                 
            end if  
        else
        
            x = 1.d0 /  Sqrt((1.d0 - bmr) * (1.d0 + bmr))
            
            alpha = (7 * r2 + b2 - 4.d0) * o9 / x
            beta = (r2 * r2 + b2 * b2 + r2 - b2 * (5.d0 + 2 * r2) + 1.d0) * o9 * x        
            gamma = bpr * x * o3 * obmr
            
            n = -4 * br * obmr * obmr
            m = 4 * br * x * x
            
            ur = 2 * r * Sqrt((1.d0 - bmr) * (1.d0 + bmr))
            ub = Sqrt((1.d0 - bmr) * (1.d0 + bmr)) * ((b + 1.d0) * (b - 1.d0) + r2) / (3 * b)
            vb = ((-1.d0 + b2)**2.d0 - 2 * (1.d0 + b2) * r2 + r2 * r2) * o3 * x / b
            
            sqomm = Sqrt(1.d0 - m)
            o = 1.d0
            ellippi = cel((sqomm), (1.d0 - n), (o), (o))
            ellipe = cel((sqomm), (o), (o), 1.d0 - m)
            ellipf = cel((sqomm), (o), (o), (o))
             
            eplusf = alpha * ellipe + beta * ellipf
            eplusf_r = ur * ellipe
            eplusf_b = ub * ellipe + vb * ellipf
            !eplusf = cel((sqomm), (o), alpha + beta, alpha * (1.d0 - m) + beta)
            !eplusf_r = cel((sqomm), (o), ur, ur * (1.d0 - m))
            !eplusf_b = cel((sqomm), (o), ub + vb, ub * (1.d0 - m) + vb)
            
            Fl = eplusf + gamma * ellippi + pisixth * (1.d0 - Sign(1.d0, bmr))
            Fl_r = eplusf_r
            Fl_b = eplusf_b
        end if
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
        Fcomplete(2) = cc * Fc_b + cl * Fl_b + cq * Fq_b
        Fcomplete(3) = cc * Fc_r + cl * Fl_r + cq * Fq_r
        Fcomplete(7) = - Fc + Fl 
        Fcomplete(8) = -2 * Fc + 2 * Fl + Fq
    else
        Fcomplete(1) = cc * Fc + cl * Fl + c2 * Fq
        Fcomplete(4) = cc * Fc_b + cl * Fl_b + c2 * Fq_b
        Fcomplete(5) = cc * Fc_r + cl * Fl_r + c2 * Fq_r
        Fcomplete(7) = - Fc + Fl 
        Fcomplete(8) = -2 * Fc + 2 * Fl + Fq
    end if
    return

end function

! evaluate the integral at one arbitrary limit along the planet or moon's boundary 
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
    real*8 :: r2, b2, tphihalf, sphihalf, cphihalf, alpha, beta, alpha_r, beta_r, gamma_r
    real*8 :: d_phi, d_r, d_b, eplusf_phi, eplusf_r, eplusf_b, lpm, lpmo
    real*8 :: ellippi, eplusf, ellipf, eplusf_m
    real*8 :: m_r, m_b, n_r, n_b, ellippim, ellippi_m, ellippi_n, ellippi_r
    real*8 :: u, v, q, p, ur, vr, qr, pr, ub, vb, qb, pb, ellippi_b, ellippe_r, kpluspi_r
    real*8 :: alpha_b, beta_b, gamma_b, nmo, mmo, mmn, sqomm, denom, cc, cl, cq
    real*8 :: ellipe
    
    r2 = r * r
    b2 = b * b
    bmr = b - r
    bpr = b + r
    br = b * r
    
    cphi = Cos(phi)
    sphi = Sin(phi)
    sphihalf = Sin(phi * 0.5)
    cphihalf = Cos(phi * 0.5)
    tphihalf = sphihalf / cphihalf
    
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
    else
        if (b == r) then
            ! This whole section (b == r) is wrong right now! 
            if (r == 0.5) then
            
                Fl = phi * o3 * 0.5 + o3 * sphihalf &
                    * (1.d0 - sphihalf * sphihalf * o3)
                Fl_phi = 0.5 * o3 * (1.d0 + 0.75 * cphihalf + 0.25 * Cos(3 * phi * 0.5))
                Fl_r = 0.d0
                Fl_b = 0.d0                
            else
                alpha = 4 * (2 * r2 - 1.d0) * o9
                alpha_r = 16 * r * o9
                
                beta = (1.d0 - 4 * r2) * o9
                beta_r = -8 * r * o9
                
                m = 4 * r2
                m_r = 8 * r
                
                nmo = n - 1.d0
                mmo = m - 1.d0
                mmn = m - n
                sqomm = Sqrt(-mmo)
                
                d = phi * o3 * 0.5 - 2 * r2 * sphi * Sqrt(1.d0 + 2 * r2 * (cphi - 1.d0)) * o9
                d_phi = 0.5 * o3 - (r2 * ((2.d0 - 4 * r2) * cphi + r2 * (1.d0 + 3*Cos(2 * phi)))) &
                        / (9 * Sqrt(1.d0 - 2 * r2 + 2 * r2 * cphi))
                d_r = (-4 * r * (1.d0 - 3 * r2 + 3 * r2 * cphi) * sphi) &
                    / (9 * Sqrt(1.d0 - 2 * r2 + 2 * r2 * cphi))
                    
                u = (mmn * nmo * (mmo * (m_r * alpha + 2 * m * alpha_r) - m_r * beta) &
                    + m * (m_r - m_r * n)) &
                    / (2 * mmo * m * mmn * nmo)
                v = (2 * beta_r * m * nmo * n - m_r * nmo * n * (alpha + beta)) &
                    / (2 * m * nmo * n)
                    
                Fl_phi = (2 * alpha - m * alpha - n * alpha + 2 * beta + m * alpha * cphi &
                    + n* alpha * cphi + n * beta * cphi + 2 * m * n * alpha * sphihalf**4) &
                    / (Sqrt(2.d0) * Sqrt(2.d0 - m + m * cphi) * (2.d0 - n + n * cphi))
                
                eplusf = el2((tphihalf), (sqomm), alpha + beta, -alpha * mmo + beta)
                eplusf_r = el2((tphihalf), (sqomm), u + v, -u * mmo + v)
                
                Fl = eplusf + d
                Fl_phi = Fl_phi + d_phi
                Fl_r = eplusf_r + d_r                
            end if            
        else if (bpr .gt. 1.d0) then
        
            y = Sqrt(br)
            
            alpha = 2 * y * (7 * r2 + b2 - 4.d0) * o9
            beta = -(3.d0 + 2*r*(b2 * b + 5 * b2 * r + 3*r*(-2.d0 + r2) + b*(-4.d0 + 7*r2))) &
                 / (18 * y)
            gamma = bpr * o3 / (2 * bmr * y)
            
            m = (1.d0 - bmr) * (1.d0 + bmr) / (4 * br)
            n = ((bmr + 1.d0) * (bmr - 1.d0)) / bmr**2.d0
            sqomm = Sqrt(1.d0 - m)

            ur = 4 * r * y
            vr = - r * (bpr + 1.d0) * (bpr - 1.d0) / y
            ub = 2 * r * (b2 + r2 - 1.d0) / (3 * y)
            vb = vr * o3
            
            if (phi_bpm .eq. 0.d0) then
            
                d = phi * o3 * 0.5 - Atan(bpr * tphihalf / bmr) * o3
                d_phi = (r2 - br * cphi) / (3 * (b2 + r2 - 2 * br * cphi))
                d_r = - b * sphi / (3 * (b2 + r2 - 2 * br * cphi))
                d_b = r * sphi / (3 * (b2 + r2 - 2 * br * cphi))
                
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
                !eplusf = cel((sqomm), (o), alpha + beta, alpha * (1.d0 - m) + beta)
                !eplusf_r = cel((sqomm), (o), ur + vr, ur * (1.d0 - m) + vr)
                !eplusf_b = cel((sqomm), (o), ub + vb, ub * (1.d0 - m) + vb)
            else                
                y = Sqrt((1.d0 - b) * (1.d0 + b) - r2 + 2 * br * cphi)
                d = o3 * (phi * 0.5 - Atan(bpr * tphihalf / bmr)) &
                    - (2 * br * o9) * sphi * y
                d_phi = o3 * o3 * r * (3 * (r - b * cphi) / (b2 + r2 - 2 * br * cphi) &
                    - 2 * b * cphi * y + 2 * b2 * r * sphi * sphi / y)
                d_r = o9 * b * (-3.d0 / (b2 + r2 - 2 * br * cphi) + 2 * (r2 - br * cphi) / y - 2 * y) * sphi
                d_b = o9 * r * (3.d0 / (b2 + r2 - 2 * br * cphi) + 2 * (b2 - br * cphi) / y - 2 * y) * sphi
                
                Fl_phi = cphihalf * (-(n * (alpha + 4 * beta)) &
                       + 4 * m * (alpha + 2*(beta + gamma)) &
                       + 4 * (m * alpha + n * beta) * cphi + n * alpha * Cos(2 * phi)) &
                       / (4 * Sqrt(m) * Sqrt(1.d0 + cphi) * Sqrt((-1.d0 + 2 * m + cphi) / m) * (2 * m - n + n * cphi))
                    
                pr = b * sphi * (3.d0 - 4 * b2 + b2 * b2 - r2 * (4.d0 + r2) + 2 * br * (4.d0 - b2 + r2) * cphi) &
                   / (9 * Sqrt(b*r) * Sqrt(2 * cphi - (b2 + r2 - 1.d0) / br) * (b2 + r2 - 2 * br * cphi))
                pb = r * sphi * (-3.d0 + 2 * b2 - b2 * b2 + 2 * r2 + r2 * r2 + 2 * br * (-2.d0 + b2 - r2) * cphi) &
                   / (9 * Sqrt(b*r) * Sqrt(2 * cphi - (b2 + r2 - 1.d0) / br) * (b2 + r2 - 2 * br * cphi))
                                        
                tans = 1.d0 / Sqrt(m / (sphihalf * sphihalf) - 1.d0)
                
                o = 1.d0
                ellippi = el3((tans), (sqomm), 1.d0 - n)
                ellipe = el2((tans), (sqomm), (o), (1.d0 - m))
                ellipf = el2((tans), (sqomm), (o), (o))
                
                eplusf = alpha * ellipe + beta * ellipf
                eplusf_r = ur * ellipe + vr * ellipf
                eplusf_b = ub * ellipe + vb * ellipf
                !eplusf = el2((tans), (sqomm), alpha + beta, alpha * (1.d0 - m) + beta)
                !eplusf_r = el2((tans), (sqomm), ur + vr, ur * (1.d0 - m) + vr)
                !eplusf_b = el2((tans), (sqomm), ub + vb, ub * (1.d0 - m) + vb)
            end if
            Fl = eplusf + gamma * ellippi + d
            Fl_phi = Fl_phi + d_phi  
            Fl_r = eplusf_r + pr + d_r
            Fl_b = eplusf_b + pb + d_b
            
        else
        
            x = Sqrt((1.d0 - bmr) * (1.d0 + bmr))
            alpha = (7 * r2 + b2 - 4.d0) * x * o9
            beta = (r2 * r2 + b2 * b2 + r2 - b2 * (5.d0 + 2 * r2) + 1.d0) / (9.d0 * x)
            gamma = bpr / (bmr * 3.d0 * x)
                        
            n = -4 * r * b / (bmr * bmr)
            m = 4 * r * b / ((1.d0 - bmr) * (1.d0 + bmr))
            
            sqomm = Sqrt(1.d0 - m)
            
            ur = 2 * r * x
            pr = b * sphi * (3.d0 - 4 * b2 + b2 * b2 - r2 * (4.d0 + r2) + 2 * br * (4.d0 - b2 + r2) * cphi) &
               / (9 * Sqrt(1.d0 - b2 - r2 + 2 * br * cphi) * (b2 + r2 - 2 * br * cphi))
                
            ub = x * ((b + 1.d0) * (b - 1.d0) + r2) / (3 * b)
            vb = (b2 * b2 + ((r - 1.d0) * (r + 1.d0))**2.d0 - 2 * b2 * (1.d0 + r2)) / (3 * b * x)
            pb = r * sphi * (-3.d0 + 2 * b2 - b2 * b2 + 2 * r2 + r2 * r2 + 2 * br * (-2.d0 + b2 - r2) * cphi) &
               / (9 * Sqrt(1.d0 - b2 - r2 + 2 * br * cphi) * (b2 + r2 - 2 * br * cphi))
            
            y = Sqrt((1.d0 - b) * (1.d0 + b) - r2 + 2 * br * cphi)
            d = o3 * (phi * 0.5 - Atan(bpr * tphihalf / bmr)) &
                - 2 * br * o9 * sphi * y
            d_r = o9 * b * sphi * (-3.d0 / (b2 + r2 - 2 * br * cphi) + 2 * (r2 - br * cphi) / y - 2 * y)
            d_b = o9 * r * (3.d0 / (b2 + r2 - 2 * br * cphi) + 2 * (b2 - br * cphi) / y - 2 * y) * sphi
            d_phi = o9 * r * (3 * (r - b * cphi) / (b2 + r2 - 2 * br * cphi) &
                - 2 * b * cphi * y + 2 * b2 * r * sphi * sphi / y)
            
            o = 1.d0
            ellippi = el3((tphihalf), (sqomm), (1.d0 - n))
            ellipe = el2((tphihalf), (sqomm), (o), (1.d0 - m))
            ellipf = el2((tphihalf), (sqomm), (o), (o))
             
            eplusf = alpha * ellipe + beta * ellipf
            eplusf_r = ur * ellipe
            eplusf_b = ub * ellipe + vb * ellipf 
            !eplusf = el2((tphihalf), (sqomm), alpha + beta, alpha * (1.d0 - m) + beta)
            !eplusf_r = el2((tphihalf), (sqomm), ur, ur * (1.d0 - m))
            !eplusf_b = el2((tphihalf), (sqomm), ub + vb, ub * (1.d0 - m) + vb)

            Fl_phi = (2 * alpha - m * alpha - n * alpha + 2 * beta - n * beta + 2 * gamma + m * alpha * cphi &
                    + n* alpha * cphi + n * beta * cphi + 2 * m * n * alpha * sphihalf**4) &
                    / (Sqrt(2.d0) * Sqrt(2.d0 - m + m * cphi) * (2.d0 - n + n * cphi))
            
            Fl = eplusf + gamma * ellippi + d
            Fl_phi = Fl_phi + d_phi
            Fl_r = eplusf_r + pr + d_r
            Fl_b = eplusf_b + pb + d_b
        end if
        
        if (b .eq. 0.d0) then
            Fl = phi * (1.d0 - (1.d0 - r2) ** (1.5)) * o3
            Fl_phi = (1.d0 - (1.d0 - r2) ** (1.5)) * o3
            Fl_r = -pi * r * Sqrt(1.d0 - r2)
        end if
        
    end if
    
2   cc = 1.d0 - c1 - 2 * c2
    cl = c1 + 2 * c2
    cq = c2
        
    Fl_bp = Fl_phi * phi_bp
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
    
    F(1) = cc * Fc + cl * Fl + cq * Fq
    F(2) = cc * Fc_bp + cl * Fl_bp + cq * Fq_bp
    F(3) = cc * Fc_rp + cl * Fl_rp + cq * Fq_rp
    F(4) = cc * Fc_bm + cl * Fl_bm + cq * Fq_bm
    F(5) = cc * Fc_rm + cl * Fl_rm + cq * Fq_rm
    F(6) = cc * Fc_bpm + cl * Fl_bpm + cq * Fq_bpm
    F(7) = - Fc + Fl 
    F(8) = -2 * Fc + 2 * Fl + Fq

    return
end function

end module phot