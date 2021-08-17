module phot
use iso_c_binding
use ellip

implicit none

real*8, parameter :: pi = 3.14159265358979323846264
real*8, parameter :: o3 = 0.33333333333333333333, o9 = 0.1111111111111111111

! inverse of full disk flux
real*8, parameter :: nc = 0.47746482927568606458024191852

contains

subroutine flux(rp, rm, bp2, bm2, bpm2, flux, j) bind(C, name="flux")

    integer (c_int), bind(C) :: j
    integer :: i
    real (c_double), bind(C) :: rp, rm
    real (c_double), bind(C) :: bp2(j), bm2(j), bpm2(j)
    real (c_double), bind(C), intent(out) :: flux(j)
    
    real*8 :: costhetapm, cosphim, cosphip, costheta
    real*8 :: thetapm, cosphim, cosphip
    real*8 :: d1, d2
    
    real*8 :: bp, bm, bpm
    bp = Sqrt(bp)
    bm = Sqrt(bm)
    bpm = Sqrt(bpm)
    
    do i=1,j,1
        
        if (bpm .gt. rp + rm) then
            if (bp .gt. rp + 1.d0) then
                if (bm .gt. rm + 1.d0) then
                    flux = 1.d0
                else
                    if (bm + rm .lt. 1.d0) then
                        flux = 1.d0 - 2 * F(pi, rm, bm) * nc
                    else
                        ! moon partially overlaps star, planet outside of star 
                    end if
                end if
            else
                if (bm .gt. rm + 1.d0) then
                    if (bp + rp .lt. 1.d0) then
                        flux = 1.d0 - 2 * F(pi, rp, bp) * nc
                    else
                        ! planet partially overlaps star, moon outside of star
                    end if
                else
                    if (bp + rp .lt. 1.d0) then
                        if (bm + rm .lt. 1.d0) then
                            flux = 1.d0 - 2 * (F(pi, rm, bm) + F(pi, rp, bp)) * nc
                        else
                            ! planet completely overlaps star, moon partially overlaps star, no mutual overlap
                        end if
                    else
                        if (bp + rp .lt. 1.d0) then
                            ! planet partially overlaps star, moon completely overlaps star, no mutual overlap
                        else
                            ! moon and planet both partially overlap star, no mutual overlap
                        end if
                    end if
                end if
            end if
        else
            if (bm .gt. rm + 1.d0) then
                if (bp + rp .lt. 1.d0) then
                    flux = 1.d0 - 2 * F(pi, rp, bp) * nc
                else
                    ! planet partially overlaps star, moon outside of star
                end if
            else
                if (bp + rp .lt. 1.d0) then
                    if (bm + rm .lt. 1.d0) then
                        if (bpm + rpm .lt. rp) then
                            flux = 1.d0 - 2 * F(pi, rp, bp) * nc
                        else
                            ! moon and planet both fully overlap star and partially overlap each other
                        end if
                    else 
                        if (bpm + rm .lt. rp) then
                            ! planet fully overlaps star, moon partially overlaps star, moon fully overlaps planet
                            ! I don't think this happens lol 
                        else
                            ! planet fully overlaps star, moon partially overlaps star, 
                            ! planet and moon partially overlap each other
                        end if
                    end if
                else
                    if (bm + rm .lt. 1.d0) then 
                        if (bpm + rm .lt. rp) then
                            ! planet partially overlaps star, moon fully overlaps star, moon fully overlaps planet
                        else
                            ! planet partially overlaps star, moon fully overlaps star, 
                            ! planet and moon partially overlap each other
                        end if
                    else
                        if (bpm + rm .lt. rp) then
                            ! planet and moon both partially overlap star, moon fully overlaps planet
                        else
                            costhetapm = (bp2 + bm2 - bpm2) / (2 * bp * bm)
                            cosphim = (bm2 + 1 - rm**2) / (2 * bm)
                            cosphip = (bp2 + 1 - rp**2) / (2 * bp)
                            thetapm = Acos(costhetapm)
                            phim = Acos(cosphim)
                            phip = Acos(cosphip)
                            if (thetapm + phim .lt. phip) then
                                ! planet and moon both partially overlap star and each other, 
                                ! but moon/star overlap is entirely overlapped by planet/star overlap
                            else if (thetapm + phip .lt. phim) then
                                ! planet and moon both partially overlap star and each other, 
                                ! but planet/star overlap is entirely overlapped by moon/star overlap
                            else
                                costheta = (bpm2 + bm2 - bp2) / (2 * bpm * bm)
                                cosphim = (bpm2 + rm**2 - rp**2) / (2 * bpm * rm)
                                cosphi1 = Cos(Acos(costheta) - Acos(cosphim))
                                cosphi2 = Cos(Acos(costheta) + Acos(cosphim))
                                d1 = rm**2 + bm2 - 2 * rm * bm * cosphi1
                                d2 = rm**2 + bm2 - 2 * rm * bm * cosphi2
                                if (d1 .gt. 1.d0) then
                                    ! planet and moon both partially overlap star and each other, 
                                    ! but the planet/moon overlap does not overlap the star
                                else if (d2 .lt. 1.d0) then
                                    ! planet and moon both partially overlap star and each other, 
                                    ! with the planet/moon overlap fully overlapping the star
                                else
                                    ! planet and moon both partially overlap star and each other, 
                                    ! with the planet/moon overlap partially overlapping the star
                                end if
                            end if
                        end if
                    end if
                end if
            end if
        end if                                    
    end do
    return
    
end do

real*8 function F(phi, r, b)

    real*8 :: phi, r, b
    real*8 :: o
    real*8 :: alpha, beta, gamma, d, s, n, m, x
    real*8 :: ellipf, ellipe, ellippi
    
    if (b == 0) then
        if (r == 0) then
            F = phi * o3
            return
        else
            F = phi * (1.d0 - (1.d0 - r * r) ** (1.5)) * o3
            return
        end if
    else if (b == r) then
        if (r == 0.5) then
            s = phi * 0.5
            F = phi * o3 * 0.5 - o3 * Sin(s) &
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
            F = alpha * ellipe + beta * ellipf + d
            return
        end if
    else if (b + r == 1.d0) then
        s = phi * 0.5
        F = phi * o3 * 0.5 - Atan((2 * r - 1.d0) / Tan(s)) * o3 &
                + Atan(2 * Sqrt(r * b) * Sin(s) / (1.d0 - 2 * r)) * o3 &
                + pi *  Sqrt(1.d0 - 4 * r * b) / (2 * r - 1.d0) * o3 * 0.5 &
                + 2.d0 * o9 * Sqrt(b * r) * (2 * r * (5 * r - 2.d0) &
                - 2 * b * r * Cos(phi) - 3.d0) * Sin(s)
        return
    else
        x = Sqrt(1.d0 - (b - r)**2.d0)
        s = phi * 0.5
        alpha = (7 * r * r + b * b - 4.d0) * x * o9
        beta = (r**4.d0 + b**4.d0 + r * r - b * b * (5.d0 + 2 * r * r) + 1.d0) / (9 * x)
        gamma = (b + r) / (b - r) / (3 * x)
        d = phi * o3 * 0.5 - Atan((b + r) / (b - r) * Tan(s)) * o3 &
                - (2 * b * r * o9) * Sin(phi) &
                * Sqrt(1.d0 - b * b - r * r + 2 * b * r * Cos(phi))
        n = - 4 * r * b / (b - r)**2
        m = 4 * r * b / (1.d0 - (r - b)**2)
        ellipf = el1(Tan(s), Sqrt(1.d0 - m))
        o = 1.d0
        ellipe = el2(Tan(s), Sqrt(1.d0 - m), o, 1.d0 - m)
        ellippi = el3(Tan(s), sqrt(1.d0 - m), 1.d0 - n)
        F = alpha * ellipe + beta * ellipf + gamma * ellippi + d
        return
    end if
        
    F = 0
    return
end function

end module phot