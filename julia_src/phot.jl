include("ellip.jl")


pihalf = 0.5 * pi; twopi = 2*pi
twopithree = twopi / 3.0
o3 = 1.0 / 3.0; o9 = 1.0 / 9.0
pithird = pi / 3.0; pisixth = pithird/2


function phis(rp::T, rm::T, bp::T, bm::T, bpm::T, cth::T, sth::T) where {T <: Real}
   #            pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta)

    #real*8 :: rp, rm, bp, bm, bpm, cth, sth
    
    # intersection angle from planet center, intersection from moon center, same 
    # values relative to bp and bm vectors respectively
    #real*8 :: pp, pm, pp1, pm1, pp2, pm2
    
    # derivatives 
    #real*8 :: pp_rp, pp_rm, pp_bpm 
    #real*8 :: pm_rp, pm_rm, pm_bpm, thetam_theta, thetam_bp, thetam_bpm
    
    # Four times the area of the triangle formed by rm, rp, and bpm
    #real*8 :: delta
    # Variables used in sorting the sides of the triangle
    #real*8 :: a, b, c, tmp
    
    # angle between bpm vector and bp vector, 
    # angle between bpm vector and bm vector
    #real*8 :: theta, thetam
    
    # for avoiding divisions
    #real*8 :: denom, obm
    
    thetam = atan(bp * sth, bpm - bp * cth)
    theta = atan(sth, cth)
    obm = 1.0 / bm
    
    # find 4 * area of triangle using modified Heron's formula 
    a = rp
    b = rm
    c = bpm
    if c > b
        tmp = c
        c = b
        b = tmp
    end
    if b > a
        tmp = b
        b = a
        a = tmp
    end
    delta = sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c)))
    denom = 1.0 / (delta * bpm * rm * rp)
    
    pm = atan(delta, (rm - rp) * (rm + rp) + bpm * bpm)   
    pm_bpm = ((rm + bpm) * (rm - bpm) - rp * rp) * denom * rm * rp
    pm_rp = 2 * rp * denom * bpm * rm * rp
    pm_rm = ((bpm - rm) * (bpm + rm) - rp * rp) * denom * rp * bpm
    
    pp = atan(delta, (rp - rm) * (rp + rm) + bpm * bpm)
    pp_bpm = ((rp - rm) * (rp + rm) - bpm * bpm) * denom * rm * rp
    pp_rp = ((bpm - rp) * (bpm + rp) - rm * rm) * denom * rm * bpm
    pp_rm = 2 * rm * denom * rm * rp * bpm
    
    thetam_bp = bpm * sth * obm * obm
    thetam_theta = ((bpm - bm) * (bpm + bm) - bp * bp) * 0.5 * obm * obm
    thetam_bpm = -bp * sth * obm * obm

    pm1 = thetam + pm
    pm2 = thetam - pm
    pp1 = theta + pp
    pp2 = theta - pp
    
    if pm1 > pi
        pm1 = pm1 - twopi
    end
    if pp1 > pi
        pp1 = pp1 - twopi
    end
    #return pm_bpm,pm_rp,pm_rm,pp_bpm,pp_rp,pp_rm,thetam_bp,thetam_theta,thetam_bpm,pm1,pm2,pp1,pp2
    return pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta
end

function kappas_p(rp::T, bp::T)  where {T <: Real}
    # , kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)

    # kp = angle to intersection from center of planet, 
    # kps = angle to intersection from center of star 
    #real*8 :: rp, bp, kp, kps
    
    # derivatives 
    #real*8 :: kp_rp, kp_bp, kps_rp, kps_bp
    
    # variables used in sorting sides of triangle
    #real*8 :: a, b, c
    
    # four times the area of the triangle with sides rp, bp, and 1
    #real*8 :: delta
    #real*8 :: denom
    
    if bp > 1.0
        a = bp
        b = 1.0
        c = rp
    else
        a = 1.0
        b = bp
        c = rp
        if rp > bp
            b = rp
            c = bp
        end
    end
    delta = sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c)))
    denom = 1.0 / (delta * bp * rp)
    
    kps = atan(delta, (1.0 - rp) * (1.0 + rp) + bp * bp)
    kps_bp = ((1.0 - bp) * (1.0 + bp) - rp * rp) * rp * denom
    kps_rp = 2 * rp * rp * bp * denom
    
    kp = atan(delta, (rp - 1.0) * (rp + 1.0) + bp * bp)
    kp_bp = ((rp + bp) * (rp - bp) - 1.0) * rp * denom
    kp_rp = ((bp + rp) * (bp - rp) - 1.0) * bp * denom
    return kp,kp_bp,kp_rp,kps,kps_bp,kps_rp
end

function kappas_m(rm::T, bp::T, bm::T, bpm::T, cth::T, sth::T) where {T <: Real}
    # , km, kms, km_rm, km_bp, km_bpm, km_theta, kms_rm, kms_bp, kms_bpm, kms_theta)

    # km = angle to interection from center of moon, 
    # kms = angle to intersection from center of planet
    #real*8 :: rm, bp, bm, bpm, cth, sth, km, kms
    
    # derivatives
    #real*8 :: km_rm, km_bp, km_bpm, km_theta, kms_rm, kms_bp, kms_bpm, kms_theta
    
    # variables used in sorting sides of triangle
    #real*8 :: a, b, c
    
    # four times the area of the triangle with sides rm, bm, and 1
    #real*8 :: delta
    
    # some useful quantities
    #real*8 :: denom, xs, xm, yp, ypm, ytheta
    
    xs = (1.0 - bm) * (1.0 + bm) - rm * rm
    xm = (rm - bm) * (rm + bm) - 1.0
    yp = bp - bpm * cth
    ypm = bpm - bp * cth
    ytheta = bp * bpm * sth
    
    if bm >= 1.0
        a = bm
        b = 1.0
        c = rm
    else
        a = 1.0
        b = bm
        c = rm
        if rm >= bm
            b = rm
            c = bm
        end
    end
    delta = sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c)))
    denom = 1.0 / (delta * bm * bm)
    
    km = atan(delta, (rm - 1.0) * (rm + 1.0) + bm * bm)
    kms = atan(delta, (1.0 - rm) * (1.0 + rm) + bm * bm)
    
    km_rm = ((bm + rm) * (bm - rm) - 1.0) / (rm * delta)
    kms_rm = 2 * rm * bm * bm * denom
    
    km_theta = ytheta * xm * denom
    kms_theta = ytheta * xs * denom
    
    km_bpm = ypm * xm * denom
    kms_bpm = ypm * xs * denom
    
    km_bp = yp * xm * denom
    kms_bp = yp * xs * denom
    return km, kms, km_rm, km_bp, km_bpm, km_theta, kms_rm, kms_bp, kms_bpm, kms_theta
end

function bm_x!(bp::T, bm::T, bpm::T, cth::T, sth::T, dbm::Array{1,T}) where {T <: Real}

    #real*8 :: bp, bm, bpm, cth, sth
    #real*8, dimension(3) :: dbm
    #real*8 :: obm 
    
    obm = 1.0 / bm
    dbm[1] = (bp - bpm * cth) * obm
    dbm[2] = (bpm - bp * cth) * obm
    dbm[3] = bp * bpm * sth * obm
    return
end 


# main loop to compute the flux at each timestep by finding the correct geometry and
# calling the integration routines 
function flux!(c1::T, c2::T, rp::T, rm::T, bp::Array{1,T}, bpm::Array{1,T}, cth::Array{1,T}, 
   sth::Array{1,T}, lc::Array{2,T}, j::Int64) where {T <: Real}

    #integer (c_int), bind(C) :: j
    #integer :: i
    #real (c_double), bind(C) :: rp, rm
    #real (c_double), bind(C), dimension(j) :: bp, cth, sth, bpm
    #real (c_double), bind(C), intent(out), dimension(8, j) :: lc
    #real*8, dimension(8) :: f0
    f0 = zeros(T,8)
    #real*8, dimension(3) :: ld
    ld = zeros(T,3)
    #real (c_double), bind(C) :: c1, c2
    #real*8 :: of0
    
    # half angles: from planet to planet/star intersection, moon to moon/star
    # intersection, spanned by planet on limb of star, spanned by moon on 
    # limb of star
    #real*8 :: kp, km, kps, kms
        
    # derivatives of above angles
    #real*8 :: kp_rp, kp_bp, kps_rp, kps_bp
    #real*8 :: km_rm, km_bp, km_bpm, km_theta
    #real*8 :: kms_rm, kms_bp, kms_bpm, kms_theta
    
    # angles to planet-moon intersection from planet center 
    #relative to bp vector (and derivatives)
    #real*8 :: pp1, pp2
    #real*8 :: pp_rp, pp_rm, pp_bpm 
    
    # angles to planet-moon intersection from moon center 
    # relative to bm vector (and derivatves)
    #real*8 :: pm1, pm2
    #real*8 :: pm_rp, pm_rm, thetam_bp, pm_bpm, thetam_theta
    # derivative of angle between bpm and bm vector with respect to bpm
    #real*8 :: thetam_bpm
    
    # used to determine cases for three body overlaps, might not be needed. 
    # Check if some of these (costheta, cosphi) can be removed when optimizing things later 
    #real*8 :: phi, phi_bpm, phi_bp, phi_bm, phi_theta, d1, d2, delta, a, b, c, tmp
    
    # For chain rule stuff
    #real*8, dimension(j) :: bm
    bm = zeros(T,j)
    #real*8, dimension(3) :: dbm, dbm0
    dbm = zeros(T,j); dbm0=zeros(T,j)
    #real*8 :: obm
    
    bm .= sqrt.((bp - bpm).^2 .+ 2 .* bp .* bpm .* (1.0 .- cth))
    dbm0 = 0.0
    
    ld[1] = 1.0 - c1 - 2 * c2
    ld[2] = c1 + 2 * c2
    ld[3] = c2
    
    # normalization factors 
    f0[1] = ld[1] * pi + ld[2] * twopithree + ld[3] * pihalf
    f0[2] = 0.0
    f0[3] = 0.0
    f0[4] = 0.0
    f0[5] = 0.0
    f0[6] = 0.0
    f0[7] = -pi + twopithree 
    f0[8] = -2 * pi + 2 * twopithree + pihalf
    
    of0 = 1.0 / f0[1]
    
    for i=1:1:j
    
        lc[:, i] = f0 * of0
    
        if bpm[i] > rp + rm
            # moon and planet don't overlap each other 
            if bp[i] < 1.0 - rp
                # planet completely inside star
                lc[:, i] .-= 2 * Fcomplete(ld, rp, bp[i], dbm0, true) * of0
            elseif (bp[i] < 1.0 + rp)
                # planet partially overlaps star
                #call kappas_p(rp, bp[i], kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                kp,kp_bp,kp_rp,kps,kps_bp,kps_rp = kappas_p(rp, bp[i])
                lc[:, i] .-= 2 * (Fstar(ld, kps, kps_rp, 0.0, kps_bp, 0.0, 0.0)
                                    .+ F(ld, kp, rp, bp[i], kp_rp, 0.0, kp_bp, 0.0, 0.0,
                                            dbm0, true, true)) * of0
            end
            if bm[i] < 1.0 - rm
                # moon completely inside star
                bm_x!(bp[i], bm[i], bpm[i], cth[i], sth[i], dbm)
                lc[:, i] .-= 2 * Fcomplete(ld, rm, bm[i], dbm, false) * of0
            elseif bm[i] < 1.0 + rm
                # moon partially overlaps star
                #call kappas_m(rm, bp[i], bm[i], bpm[i], cth[i], sth[i], km, kms,
                #              km_rm, km_bp, km_bpm, km_theta,
                #              kms_rm, kms_bp, kms_bpm, kms_theta)
                km, kms, km_rm, km_bp, km_bpm, km_theta, kms_rm, kms_bp, kms_bpm, kms_theta = kappas_m(rm, bp[i], bm[i], bpm[i], cth[i], sth[i])
                #call bm_x(bp[i], bm[i], bpm[i], cth[i], sth[i], dbm)
                bm_x!(bp[i], bm[i], bpm[i], cth[i], sth[i], dbm)
                lc[:, i] .-= 2 * (Fstar(ld, kms, 0.0, kms_rm, kms_bp, kms_bpm, kms_theta)
                                    .+ F(ld, km, rm, bm[i], 0.0, km_rm, km_bp, km_bpm, km_theta,
                                            dbm, false, true)) * of0
            end
        else
            # moon and planet do overlap each other 
            if bp[i] > rp + 1.0
                if bm[i] > rm + 1.0
                    # neither moon nor planet overlap star
                    lc[:, i] .= f0 .* of0
                else
                    # moon partially overlaps star, planet does not overlap star
                    #call bm_x(bp[i], bm[i], bpm[i], cth[i], sth[i], dbm)
                    bm_x!(bp[i], bm[i], bpm[i], cth[i], sth[i], dbm)
                    #call kappas_m(rm, bp[i], bm[i], bpm[i], cth[i], sth[i], km, kms, &
                    #                  km_rm, km_bp, km_bpm, km_theta, &
                    #                  kms_rm, kms_bp, kms_bpm, kms_theta)
                    km, kms, km_rm, km_bp, km_bpm, km_theta, kms_rm, kms_bp, kms_bpm, kms_theta = kappas_m(rm, bp[i], bm[i], bpm[i], cth[i], sth[i])
                    lc[:, i] .= 2 * (Fstar(ld, pi - kms, 0.0, -kms_rm, -kms_bp, -kms_bpm, -kms_theta)
                                    .- F(ld, km, rm, bm[i], 0.0, km_rm, km_bp, km_bpm, km_theta,
                                            dbm, false, true)) * of0
                end
            else
                if bm[i] > rm + 1.0
                    # planet partially overlaps star, moon is outside of star
                    #call kappas_p(rp, bp[i], kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                    kp,kp_bp,kp_rp,kps,kps_bp,kps_rp = kappas_p(rp, bp[i])
                    lc[:, i] .= 2 * (Fstar(ld, pi - kps, -kps_rp, 0.0, -kps_bp, 0.0, 0.0)
                                    .- F(ld, kp, rp, bp[i], kp_rp, 0.0, kp_bp, 0.0, 0.0,
                                        dbm0, true, true)) * of0
                else
                    if bp[i] + rp <= 1.0
                        if bm[i] + rm <= 1.0
                            if bpm[i] + rm <= rp
                                # moon and planet both overlap star, moon fully overlapped by planet
                                lc[:, i] .= (f0 - 2 * Fcomplete(ld, rp, bp[i], dbm0, true)) * of0
                            else
                                # Case E
                                # moon and planet both overlap star, moon and planet partially overlap each other 
                                #call bm_x(bp[i], bm[i], bpm[i], cth[i], sth[i], dbm)
                                bm_x!(bp[i], bm[i], bpm[i], cth[i], sth[i], dbm)
                                #call phis(rp, rm, bp[i], bm[i], bpm[i], cth[i], sth[i], pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, &
                                #          pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta)
                                pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta =
                                   phis(rp, rm, bp[i], bm[i], bpm[i], cth[i], sth[i])
                                lc[:, i] .= (f0 .- Arc(ld, pp1, pp2, rp, bp[i],
                                                     pp_rp, pp_rm, 0.0, pp_bpm, 1.0,
                                                     -pp_rp, -pp_rm, 0.0, -pp_bpm, 1.0,
                                                     dbm0, true, false, false)
                                                .- Arc(ld, pm1, pm2, rm, bm[i],
                                                     pm_rp, pm_rm, thetam_bp, pm_bpm + thetam_bpm, thetam_theta,
                                                     -pm_rp, -pm_rm, thetam_bp, -pm_bpm + thetam_bpm, thetam_theta,
                                                     dbm, false, false, false)) * of0
                            end
                        else
                            # planet fully overlaps star, moon partially overlaps star, both overlap each other 
                            #call bm_x(bp[i], bm[i], bpm[i], cth[i], sth[i], dbm)
                            bm_x!(bp[i], bm[i], bpm[i], cth[i], sth[i], dbm)
                            #call phis(rp, rm, bp[i], bm[i], bpm[i], cth[i], sth[i], pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm,
                            #          pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta)
                            pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta =
                                   phis(rp, rm, bp[i], bm[i], bpm[i], cth[i], sth[i])
                            #call kappas_m(rm, bp[i], bm[i], bpm[i], cth[i], sth[i], km, kms,
                            #          km_rm, km_bp, km_bpm, km_theta, &
                            #          kms_rm, kms_bp, kms_bpm, kms_theta)
                            km, kms, km_rm, km_bp, km_bpm, km_theta, kms_rm, kms_bp, kms_bpm, kms_theta = kappas_m(rm, bp[i], bm[i], bpm[i], cth[i], sth[i])
                            lc[:, i] .= (2 * Fstar(ld, pi - kms, 0.0, -kms_rm, -kms_bp, -kms_bpm, -kms_theta)
                                        .- Arc(ld, -km, pm2, rm, bm[i], 
                                              0.0, -km_rm, -km_bp, -km_bpm, -km_theta,
                                              -pm_rp, -pm_rm, thetam_bp, -pm_bpm + thetam_bpm, thetam_theta,
                                              dbm, false, true, false)
                                        .- Arc(ld, pm1, km, rm, bm[i],
                                              pm_rp, pm_rm, thetam_bp, pm_bpm + thetam_bpm, thetam_theta,
                                              0.0, km_rm, km_bp, km_bpm, km_theta,
                                              dbm, false, false, true)
                                        .- Arc(ld, pp1, pp2, rp, bp[i],
                                              pp_rp, pp_rm, 0.0, pp_bpm, 1.0, 
                                              -pp_rp, -pp_rm, 0.0, -pp_bpm, 1.0, 
                                              dbm0, true, false, false)) * of0
                        end
                    else
                        if bm[i] + rm <= 1.0
                            if bpm[i] + rm <= rp
                                # planet partially overlaps star, moon fully overlaps star but is completely overlapped by planet 
                                #call kappas_p(rp, bp[i], kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                                kp,kp_bp,kp_rp,kps,kps_bp,kps_rp = kappas_p(rp, bp[i])
                                lc[:, i] = 2 * (Fstar(ld, pi - kps, -kps_rp, 0.0, -kps_bp, 0.0, 0.0)
                                                .- F(ld, kp, rp, bp[i], kp_rp, 0.0, kp_bp, 0.0, 0.0,
                                                    dbm0, true, true)) * of0
                            else
                                # planet partially overlaps star, moon fully overlaps star and only partially overlaps planet
                                #call bm_x(bp[i], bm[i], bpm[i], cth[i], sth[i], dbm)
                                bm_x!(bp[i], bm[i], bpm[i], cth[i], sth[i], dbm)
                                #call phis(rp, rm, bp[i], bm[i], bpm[i], cth[i], sth[i], pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, &
                                #      pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta)
                                pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta =
                                   phis(rp, rm, bp[i], bm[i], bpm[i], cth[i], sth[i])
                                #call kappas_p(rp, bp[i], kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                                kp,kp_bp,kp_rp,kps,kps_bp,kps_rp = kappas_p(rp, bp[i])
                                lc[:, i] .= (2 * Fstar(ld, pi - kps, -kps_rp, 0.0, -kps_bp, 0.0, 0.0)
                                            .- Arc(ld, -kp, pp2, rp, bp[i],
                                                  -kp_rp, 0.0, -kp_bp, 0.0, 0.0,
                                                  -pp_rp, -pp_rm, 0.0, -pp_bpm, 1.0,
                                                  dbm0, true, true, false)
                                            .- Arc(ld, pp1, kp, rp, bp[i],
                                                  pp_rp, pp_rm, 0.0, pp_bpm, 1.0,
                                                  kp_rp, 0.0, kp_bp, 0.0, 0.0,
                                                  dbm0, true, false, true)
                                            .- Arc(ld, pm1, pm2, rm, bm[i],
                                                  pm_rp, pm_rm, thetam_bp, pm_bpm + thetam_bpm, thetam_theta,
                                                  -pm_rp, -pm_rm, thetam_bp, -pm_bpm + thetam_bpm, thetam_theta,
                                                  dbm, false, false, false)) * of0
                            end
                        else
                            if bpm[i] + rm <= rp
                                # planet and moon both partially overlap star but moon is fully overlapped by the planet
                                #call kappas_p(rp, bp[i], kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                                kp,kp_bp,kp_rp,kps,kps_bp,kps_rp = kappas_p(rp, bp[i])
                                lc[:, i] .= 2 * (Fstar(ld, pi - kps, -kps_rp, 0.0, -kps_bp, 0.0, 0.0)
                                                - F(ld, kp, rp, bp[i], kp_rp, 0.0, kp_bp, 0.0, 0.0,
                                                    dbm0, true, true)) * of0
                            else
                                #call kappas_p(rp, bp[i], kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                                kp,kp_bp,kp_rp,kps,kps_bp,kps_rp = kappas_p(rp, bp[i])
                                #call kappas_m(rm, bp[i], bm[i], bpm[i], cth[i], sth[i], km, kms, &
                                #      km_rm, km_bp, km_bpm, km_theta, &
                                #      kms_rm, kms_bp, kms_bpm, kms_theta)
                                km, kms, km_rm, km_bp, km_bpm, km_theta, kms_rm, kms_bp, kms_bpm, kms_theta = kappas_m(rm, bp[i], bm[i], bpm[i], cth[i], sth[i])
                                
                                a = bm[i]
                                b = bp[i]
                                c = bpm[i]
                                if b > a
                                    tmp = b
                                    b = a
                                    a = tmp
                                end
                                if c > b
                                    tmp = c
                                    c = b
                                    b = tmp
                                end
                                if b > a
                                    tmp = b
                                    b = a
                                    a = tmp
                                end
                                delta = sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c)))
                                phi = atan(delta, (bm[i] - bpm[i]) * (bm[i] + bpm[i]) + bp[i] * bp[i])                                
                                
                                #call phis(rp, rm, bp[i], bm[i], bpm[i], cth[i], sth[i], pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, &
                                #          pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta)
                                pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta =
                                   phis(rp, rm, bp[i], bm[i], bpm[i], cth[i], sth[i])
                                
                                if phi + kms <= kps
                                        if pp2 > kp
                                            # planet and moon both partially overlap the star and each other but the 
                                            # moon-star overlap is contained within the planet-star overlap
                                            lc[:, i] = 2 * (Fstar(ld, pi - kps, -kps_rp, 0.0, -kps_bp, 0.0, 0.0)
                                                            .- F(ld, kp, rp, bp[i], kp_rp, 0.0, kp_bp, 0.0, 0.0,
                                                                dbm0, true, true)) * of0
                                        else
                                            # planet and moon both partially overlap star and each other but the 
                                            # moon-star intersections are overlapped by the planet
                                            #call bm_x(bp[i], bm[i], bpm[i], cth[i], sth[i], dbm)
                                            bm_x!(bp[i], bm[i], bpm[i], cth[i], sth[i], dbm)
                                            lc[:, i] = (2 * Fstar(ld, pi - kps, -kps_rp, 0.0, -kps_bp, 0.0, 0.0)
                                                        .- Arc(ld, -kp, pp2, rp, bp[i],
                                                              -kp_rp, 0.0, -kp_bp, 0.0, 0.0,
                                                              -pp_rp, -pp_rm, 0.0, -pp_bpm, 1.0,
                                                              dbm0, true, true, false)
                                                        .- Arc(ld, pp1, kp, rp, bp[i],
                                                              pp_rp, pp_rm, 0.0, pp_bpm, 1.0,
                                                              kp_rp, 0.0, kp_bp, 0.0, 0.0,
                                                              dbm0, true, false, true)
                                                        .- Arc(ld, pm1, pm2, rm, bm[i],
                                                              pm_rp, pm_rm, thetam_bp, pm_bpm + thetam_bpm, thetam_theta,
                                                              -pm_rp, -pm_rm, thetam_bp, -pm_bpm + thetam_bpm, thetam_theta,
                                                              dbm, false, false, false)) * of0
                                        end
                                elseif phi + kps <= kms
                                    #call bm_x(bp[i], bm[i], bpm[i], cth[i], sth[i], dbm)
                                    bm_x!(bp[i], bm[i], bpm[i], cth[i], sth[i], dbm)
                                    if ((bp[i] - rp) <=  (bm[i] - rm))
                                        # planet and moon both partially overlap the star and each other but the 
                                        # planet-star intersections are overlapped by the moon
                                        lc[:, i] .= (2 * Fstar(ld, pi - kms, 0.0, -kms_rm, -kms_bp, -kms_bpm, -kms_theta)
                                                        .- Arc(ld, -km, pm2, rm, bm[i],
                                                              0.0, -km_rm, -km_bp, -km_bpm, -km_theta,
                                                              -pm_rp, -pm_rm, thetam_bp, -pm_bpm + thetam_bpm, thetam_theta,
                                                              dbm, false, true, false)
                                                        .- Arc(ld, pm1, km, rm, bm[i],
                                                              pm_rp, pm_rm, thetam_bp, pm_bpm + thetam_bpm, thetam_theta,
                                                              0.0, km_rm, km_bp, km_bpm, km_theta,
                                                              dbm, false, false, true)
                                                        .- Arc(ld, pp1, pp2, rp, bp[i],
                                                              pp_rp, pp_rm, 0.0, pp_bpm, 1.0,
                                                              -pp_rp, -pp_rm, 0.0, -pp_bpm, 1.0,
                                                              dbm0, true, false, false)) * of0
                                    else
                                        # planet and moon both partially overlap the star and each other but 
                                        # the planet-star overlap is  entirely within the moon-star overlap
                                        lc[:, i] .= 2 * (Fstar(ld, pi - kms, 0.0, -kms_rm, -kms_bp, -kms_bpm, -kms_theta)
                                                        .- F(ld, km, rm, bm[i], 0.0, km_rm, km_bp, km_bpm, km_theta,
                                                            dbm, false, true)) * of0
                                    end
                                else
                                    # bookmark
                                    d1 = rm * rm + bm[i] * bm[i] - 2 * rm * bm[i] * cos(pm2)
                                    d2 = rm * rm + bm[i] * bm[i] - 2 * rm * bm[i] * cos(pm1)
                                    #call bm_x(bp[i], bm[i], bpm[i], cth[i], sth[i], dbm)
                                    bm_x!(bp[i], bm[i], bpm[i], cth[i], sth[i], dbm)
                                    if ((d1 > 1.0) && (d2 > 1.0))
                                        # planet and moon both partially overlap star and each other, 
                                        # but the planet/moon overlap does not overlap the star
                                        lc[:, i] .= 2 * (Fstar(ld, pi - (kps + kms), -kps_rp, -kms_rm,
                                                              -(kps_bp + kms_bp), -kms_bpm, -kms_theta)
                                                        .- F(ld, kp, rp, bp[i], kp_rp, 0.0, kp_bp, 0.0, 0.0,
                                                            dbm0, true, true)
                                                        .- F(ld, km, rm, bm[i], 0.0, km_rm, km_bp, km_bpm, km_theta,
                                                            dbm, false, true)) * of0
                                    elseif ((d1 <= 1.0) && (d2 <= 1.0))
                                        # planet and moon both partially overlap star and each other, 
                                        # with the planet/moon overlap fully overlapping the star
                                        lc[:, i] .= (2 * Fstar(ld, pi - (kps + kms), -kps_rp, -kms_rm,
                                                              -(kps_bp + kms_bp), -kms_bpm, -kms_theta)
                                                        .- Arc(ld, -km, -pm1, rm, bm[i],
                                                              0.0, -km_rm, -km_bp, -km_bpm, -km_theta,
                                                              -pm_rp, -pm_rm, -thetam_bp, -pm_bpm - thetam_bpm, -thetam_theta,
                                                              dbm, false, true, false)
                                                        .- Arc(ld, -pm2, km, rm, bm[i],
                                                              pm_rp, pm_rm, -thetam_bp, pm_bpm - thetam_bpm, -thetam_theta,
                                                              0.0, km_rm, km_bp, km_bpm, km_theta,
                                                              dbm, false, false, true)
                                                        .- Arc(ld, pp1, kp, rp, bp[i],
                                                              pp_rp, pp_rm, 0.0, pp_bpm, 1.0,
                                                              kp_rp, 0.0, kp_bp, 0.0, 0.0,
                                                              dbm0, true, false, true)
                                                        .- Arc(ld, -kp, pp2, rp, bp[i],
                                                              -kp_rp, 0.0, -kp_bp, 0.0, 0.0,
                                                              -pp_rp, -pp_rm, 0.0, -pp_bpm, 1.0,
                                                              dbm0, true, true, false)) * of0
                                    else
                                        # planet and moon both partially overlap star and each other, 
                                        # with the planet/moon overlap partially overlapping the star
                                        
                                        # there might be a mistake somewhere in here... 
                                        phi_bpm = bp[i] * sth[i] / (bm[i] * bm[i])
                                        phi_theta = bpm[i] * (bp[i] * cth[i] - bpm[i]) / (bm[i] * bm[i])
                                        phi_bp = - bpm[i] * sth[i] / (bm[i] * bm[i])
                                        
                                        lc[:, i] .= (2 * Fstar(ld, pi - 0.5 * (kps + kms + phi),
                                                              -0.5 * kps_rp, -0.5 * kms_rm,
                                                              -0.5 * (kps_bp + kms_bp + phi_bp), -0.5 * (kms_bpm + phi_bpm),
                                                              -0.5 * (kms_theta + phi_theta))
                                                        .- Arc(ld, -pm2, km, rm, bm[i],
                                                              pm_rp, pm_rm, -thetam_bp, pm_bpm - thetam_bpm, -thetam_theta,
                                                              0.0, km_rm, km_bp, km_bpm, km_theta,
                                                              dbm, false, false, true)
                                                        .- Arc(ld, -kp, pp2, rp, bp[i],
                                                              -kp_rp, 0.0, -kp_bp, 0.0, 0.0,
                                                              -pp_rp, -pp_rm, 0.0, -pp_bpm, 1.0,
                                                              dbm0, true, true, false)) * of0
                                    end 
                                end
                            end
                        end
                    end
                end
            end
        end
        lc[:, i] .-=  f0 .* of0
    end
    return
end

# work out the right sign and order of the integration and call the integration routine 
# to integrate along an arbitrary arc of the planet or moon 
function Arc(ld::Array{1,T}, phi1::T, phi2::T, r::T, b::T, phi1_rp::T, 
                phi1_rm::T, phi1_bp::T, phi1_bpm::T, phi1_theta::T, phi2_rp::T, phi2_rm::T, 
                phi2_bp::T, phi2_bpm::T, phi2_theta::T, dbm::Array{1,T}, 
                pflag::Bool, limbflag1::Bool, limbflag2::Bool)  where {T <: Real}
                    
    #real*8, dimension(8) :: Arc
    Arc = zeros(T,8)
    #logical :: pflag, limbflag1, limbflag2
    #real*8 :: phi1, phi2, r, b
    #real*8, dimension(3) :: ld
    #real*8, dimension(3) :: dbm
    #real*8 :: phi1_rp, phi1_rm, phi1_bp, phi1_bpm, phi1_theta
    #real*8 :: phi2_rp, phi2_rm, phi2_bp, phi2_bpm, phi2_theta
    #real*8 :: const, lin, quad
        
    if (phi1 < 0)
        if (phi2 > 0)
            Arc .= (F(ld, phi2, r, b, phi2_rp, phi2_rm, phi2_bp, phi2_bpm, phi2_theta,
                    dbm, pflag, limbflag2)
                .+ F(ld, -phi1, r, b, -phi1_rp, -phi1_rm, -phi1_bp, -phi1_bpm, -phi1_theta,
                    dbm, pflag, limbflag1))
            return Arc
        else
            if (phi2 < phi1)
                Arc .= (2 * Fcomplete(ld, r, b, dbm, pflag)
                    .+ F(ld, -phi1, r, b, -phi1_rp, -phi1_rm, -phi1_bp, -phi1_bpm, -phi1_theta,
                        dbm, pflag, limbflag1)
                    .- F(ld, -phi2, r, b, -phi2_rp, -phi2_rm, -phi2_bp, -phi2_bpm, -phi2_theta,
                        dbm, pflag, limbflag2))
                return Arc
            else
                Arc .= (- F(ld, -phi2, r, b, -phi2_rp, -phi2_rm, -phi2_bp, -phi2_bpm, -phi2_theta,
                          dbm, pflag, limbflag2)
                      .+ F(ld, -phi1, r, b, -phi1_rp, -phi1_rm, -phi1_bp, -phi1_bpm, -phi1_theta,
                          dbm, pflag, limbflag1))
                return Arc
            end
        end
    else
        if (phi2 < 0)
            Arc .= (2 * Fcomplete(ld, r, b, dbm, pflag)
                .- F(ld, phi1, r, b, phi1_rp, phi1_rm, phi1_bp, phi1_bpm, phi1_theta,
                    dbm, pflag, limbflag1)
                .- F(ld, -phi2, r, b, -phi2_rp, -phi2_rm, -phi2_bp, -phi2_bpm, -phi2_theta,
                    dbm, pflag, limbflag2))
            return Arc
        else
            if (phi2 < phi1)
                Arc .= (2 * Fcomplete(ld, r, b, dbm, pflag)
                    .+ F(ld, phi2, r, b, phi2_rp, phi2_rm, phi2_bp, phi2_bpm, phi2_theta,
                        dbm, pflag, limbflag2)
                    .- F(ld, phi1, r, b, phi1_rp, phi1_rm, phi1_bp, phi1_bpm, phi1_theta,
                        dbm, pflag, limbflag1))
                return Arc
            else
                Arc .= (F(ld, phi2, r, b, phi2_rp, phi2_rm, phi2_bp, phi2_bpm, phi2_theta,
                        dbm, pflag, limbflag2)
                    .- F(ld, phi1, r, b, phi1_rp, phi1_rm, phi1_bp, phi1_bpm, phi1_theta,
                        dbm, pflag, limbflag1))
                return Arc
            end
        end
    end
    
    return Arc
    
end

# integrate along the limb of the star
function Fstar(ld::Array{1,T}, phi::T, phi_rp::T, phi_rm::T, phi_bp::T, phi_bpm::T, phi_theta::T) where {T <: Real}

    #real*8, dimension(8) :: Fstar
    Fstar = zeros(T,8)
    #real*8, dimension(3) :: F_, F_rp, F_rm, F_bp, F_bpm, F_theta
    F_=zeros(T,3); F_rp=zeros(3); F_rm=zeros(3); F_bp=zeros(3); F_bpm=zeros(3); F_theta=zeros(3)

    #real*8 :: phi
    #real*8 :: Fc_phi, Fq_phi, Fl_phi
    #real*8 :: phi_bp, phi_rp, phi_bm, phi_rm, phi_bpm, phi_theta
    #real*8, dimension(3) :: ld
    #ld = zeros(T,3)
        
    F_[1] = 0.5 * phi
    Fc_phi = 0.5
    F_rp[1] = Fc_phi * phi_rp
    F_rm[1] = Fc_phi * phi_rm
    F_bp[1] = Fc_phi * phi_bp
    F_bpm[1] = Fc_phi * phi_bpm
    F_theta[1] = Fc_phi * phi_theta
    
    F_[2] = phi * o3
    Fl_phi = o3
    F_rp[2] = Fl_phi * phi_rp
    F_rm[2] = Fl_phi * phi_rm
    F_bp[2] = Fl_phi * phi_bp
    F_bpm[2] = Fl_phi * phi_bpm
    F_theta[2] = Fl_phi * phi_theta

    #F_[3] = 0.25 * (phi - o3 * sin(phi) * cos(3 * phi))
    #Fq_phi = 0.5 * o3 * cos(phi)^2 * (5.0 - 4 * cos(2 * phi))
    F_[3] = 0.25 * phi
    Fq_phi = 0.25
    F_rp[3] = Fq_phi * phi_rp
    F_rm[3] = Fq_phi * phi_rm
    F_bp[3] = Fq_phi * phi_bp
    F_bpm[3] = Fq_phi * phi_bpm
    F_theta[3] = Fq_phi * phi_theta
    
    Fstar[1] = sum(ld .* F_)
    Fstar[2] = sum(ld .* F_rp)
    Fstar[3] = sum(ld .* F_rm)
    Fstar[4] = sum(ld .* F_bp)
    Fstar[5] = sum(ld .* F_bpm)
    Fstar[6] = sum(ld .* F_theta)
    Fstar[7] = - F_[1] + F_[2]
    Fstar[8] = -2 * F_[1] + 2 * F_[2] + F_[3]
    return Fstar
end

# integrate around the entire planet/moon 
function Fcomplete(ld::Array{1,T}, r::T, b::T, dbm::Array{1,T}, pflag::Bool) where {T <: Real}

    #real*8, dimension(8) :: Fcomplete
    Fcomplete = zeros(T,8)
    #real*8, dimension(3) :: F_, F_r, F_b
    F_ = zeros(T,3); F_r = zeros(T,3); F_b = zeros(T,3)
    #real*8 :: sumF_b
    
    # Are we integrating along the edge of the moon or the planet? 
    #logical :: pflag
    
    # Limb darkening params
    #real*8, dimension(3) :: ld
    
    # self explanatory
    #real*8 :: r, b
    
    # derivatives of input parameters 
    #real*8, dimension(3) :: dbm
    
    # convenient parameters
    #real*8 :: r2, b2, br, bmr, bpr, obmr
    #real*8 :: x, y, ox, ome, o
    
    # components of flux and their derivatives
    #real*8 :: Fc, Fc_r, Fc_b
    #real*8 :: Fq, Fq_r, Fq_b
    #real*8 :: Fl, Fl_r, Fl_b
    
    # For the integral
    #real*8 :: alpha, beta, gamma, d, n, m
    #real*8 :: d_r, d_b, sgn
    #real*8 :: ur, vr, ub, vb, sqomm
    
    # Elliptic integrals
    #real*8 :: ellippi, ellipe, ellipf
    #real*8 :: eplusf, eplusf_r, eplusf_b
    
    r2 = r * r
    b2 = b * b
    bmr = b - r
    obmr = 1.0 / bmr
    bpr = b + r
    br = b * r
    
    F_[1] = r2 * pihalf  
    F_r[1] = r * pi
    F_b[1] = 0.0
    
    F_[3] = pihalf * r2 * (b2 + 0.5 * r2) 
    F_r[3] = pi * r * (b2 + r2)
    F_b[3] = pi * r2 * b
    
    if ld[2] == 0.0
        Fl = 0.0
        Fl_r = 0.0
        Fl_b = 0.0        
    else
        x = sqrt((1.0 - bmr) * (1.0 + bmr))
        ox = 1.0 / x
            
        alpha = (7 * r2 + b2 - 4.0) * o9 * x
        beta = (r2 * r2 + b2 * b2 + r2 - b2 * (5.0 + 2 * r2) + 1.0) * o9 * ox        
            
        n = -4 * br * obmr * obmr
        m = 4 * br * ox * ox
            
        ur = 2 * r * x
        ub = x * ((b + 1.0) * (b - 1.0) + r2) / (3 * b)
        vb = ((-1.0 + b2)^2 - 2 * (1.0 + b2) * r2 + r2 * r2) * o3 * ox / b
            
        sqomm = sqrt(1.0 - m)
        o = 1.0
            
        if b == r
            ellippi = 0.0
            sgn = 0.0
            gamma = 0.0
        else
            ellippi = cel((sqomm), (1.0 - n), (o), (o))
            sgn = sign(bmr)
            gamma = bpr * ox * o3 * obmr
        end
            
        ellipe = cel((sqomm), (o), (o), 1.0 - m)
        ellipf = cel((sqomm), (o), (o), (o))
             
        eplusf = alpha * ellipe + beta * ellipf
        eplusf_r = ur * ellipe
        eplusf_b = ub * ellipe + vb * ellipf
            
        F_[2] = eplusf + gamma * ellippi + pisixth * (1.0 - sgn)
        F_r[2] = eplusf_r
        F_b[2] = eplusf_b
    end

    Fcomplete .= 0.0
    sumF_b = sum(ld .* F_b)
    
    if pflag
        Fcomplete[1] = sum(ld .* F_)
        Fcomplete[2] = sum(ld .* F_r)
        Fcomplete[4] = sum(ld .* F_b)
        Fcomplete[7] = - F_[1] + F_[2]
        Fcomplete[8] = -2 * F_[1] + 2 * F_[2] + F_[3]
    else
        Fcomplete[1] = sum(ld .* F_)
        Fcomplete[3] = sum(ld .* F_r)
        Fcomplete[4] = sumF_b .* dbm[1]
        Fcomplete[5] = sumF_b .* dbm[2]
        Fcomplete[6] = sumF_b .* dbm[3]
        Fcomplete[7] = - F_[1] + F_[2]
        Fcomplete[8] = -2 * F_[1] + 2 * F_[2] + F_[3]
    end
    return Fcomplete
end

# evaluate the integral at one arbitrary limit along the planet or moon's boundary 
function F(ld::Array{1,T}, phi::T, r::T, b::T, phi_rp::T, phi_rm::T, 
     phi_bp::T, phi_bpm::T, phi_theta::T, dbm::Array{1,T}, pflag::Bool, limbflag::Bool) where {T <: Real}

    #real*8, dimension(8) :: F
    F = zeros(T,8)
    #real*8, dimension(3) :: F_, F_rp, F_rm, F_bp, F_bpm, F_theta
    F_ = zeros(T,3); F_rp= zeros(T,3); F_rm= zeros(T,3); F_bp= zeros(T,3); F_bpm= zeros(T,3); F_theta= zeros(T,3); 
    
    # Are we integrating along the edge of the moon or the planet? 
    #logical :: pflag, limbflag
    
    # Limb darkening params
    #real*8, dimension(3) :: ld
    
    # self explanatory
    #real*8 :: phi, r, b
    
    # derivatives of input parameters 
    #real*8 :: phi_bp, phi_rp, phi_bm, phi_rm, phi_bpm, phi_theta
    #real*8, dimension(3) :: dbm
    
    # convenient parameters
    #real*8 :: sphi, cphi, tphihalf, sphihalf, cphihalf
    #real*8 :: r2, b2, br, bmr, bpr, obmr
    #real*8 :: x, y, z, ox, oy, oz, ome, tans, o
    
    # components of flux and their derivatives
    #real*8 :: Fc, Fc_phi, Fc_r, Fc_b
    #real*8 :: Fq, Fq_phi, Fq_r, Fq_b
    #real*8 :: Fl, Fl_phi, Fl_r, Fl_b
    
    # For the integral
    #real*8 :: alpha, beta, gamma, d, n, m
    #real*8 :: d_phi, d_r, d_b
    #real*8 :: ur, vr, ub, vb, pr, pb, sqomm
    
    # Elliptic integrals
    #real*8 :: ellippi, ellipe, ellipf
    #real*8 :: eplusf, eplusf_r, eplusf_b
    
    r2 = r * r
    b2 = b * b
    bmr = b - r
    obmr = 1.0 / bmr
    bpr = b + r
    br = b * r
    
    cphi = cos(phi)
    sphi = sin(phi)
    sphihalf = sin(phi * 0.5)
    cphihalf = cos(phi * 0.5)
    tphihalf = sphihalf / cphihalf
    
    if phi == pi
        F .= Fcomplete(ld, r, b, dbm, pflag)
        return F
    end
        
    F_[1] = 0.5 * (r2 * phi - br * sphi)
    Fc_phi = 0.5 * (r2 - br * cphi)
    Fc_b = -0.5 * r * sphi
    Fc_r = 0.5 * (2 * r * phi - b * sphi)
    
    F_rp[1] = Fc_phi * phi_rp
    F_rm[1] = Fc_phi * phi_rm
    F_bp[1] = Fc_phi * phi_bp
    F_bpm[1] = Fc_phi * phi_bpm
    F_theta[1] = Fc_phi * phi_theta
    
    #F_[3] = -0.25 * 0.25 * o3 * (r * (4 * b * (2 * b2 + 9 * r2) * sphi 
    #   - 4 * r * (3 * (2 * b2 + r2) * phi 
    #   + br * sin(3 * phi)) + r2 * r * sin(4 * phi)))
    F_[3] = 0.25 * (r * (r * (2 * b2 + r2) * phi 
           + b * (-b2 - 3 * r2 + br * cphi) * sphi))
    #Fq_phi = - o3 * 0.25 * (r * (b * (2 * b2 + 9 * r2) * cphi
    #    - 3 * r * (2 * b2 + r2 + br * cos(3 * phi)) + r2 * r * cos(4 * phi)))
    Fq_phi = 0.25 * (r * (2 * b2 * r + r2 * r + b * (-((b2 + 3 * r2) * cphi)
            + br * Cos(2 * phi))))
    #Fq_r = ( - (b * (2 * b2 + 27 * r2) * sphi)
    #     + r * (12 * (b2 + r2) * phi + 3 * br * sin(3 * phi)
    #     - r2 * sin(4 * phi))) * 0.25 * o3
    Fq_r = r * (b2 + r2) * phi - 0.25 * (b * (b2 + 9 * r2 - 2 * br * cphi) * sphi)

    #Fq_b = (b * (r2 * phi - br * sphi * 0.5)
    #    + r2 * 0.25 * r * (sin(3 * phi) * o3 - 3 * sphi))
    Fq_b = 0.25 * (r * (-3 * (b2 + r2) * sphi + br * (4 * phi + Sin(2 * phi))))
            
    F_rp[3] = Fq_phi * phi_rp
    F_rm[3] = Fq_phi * phi_rm
    F_bp[3] = Fq_phi * phi_bp
    F_bpm[3] = Fq_phi * phi_bpm
    F_theta[3] = Fq_phi * phi_theta
        
    if ld[2] == 0.0
        Fl = 0.0
        Fl_phi = 0.0
        Fl_r = 0.0
        Fl_b = 0.0        
    else          
        if bpr > 1.0
        
            y = sqrt(br)
            oy = 1.0 / y
            ox = 1.0 / (b2 + r2 - 2 * br * cphi)
            
            alpha = 2 * y * (7 * r2 + b2 - 4.0) * o9
            beta = -(3.0 + 2*r*(b2 * b + 5 * b2 * r + 3*r*(-2.0 + r2) + b*(-4.0 + 7*r2))) * o9 * 0.5 * oy
            gamma = bpr * o3 * 0.5 * oy * obmr
            
            m = (1.0 - bmr) * (1.0 + bmr) * 0.25 * oy * oy
            n = ((bmr + 1.0) * (bmr - 1.0)) * obmr * obmr
            sqomm = sqrt(1.0 - m)

            ur = 4 * r * y
            vr = - r * (bpr + 1.0) * (bpr - 1.0) * oy
            ub = 2 * r * (b2 + r2 - 1.0) * o3 * oy
            vb = vr * o3
            
            if limbflag
            
                d = phi * o3 * 0.5 - atan(bpr * tphihalf * obmr) * o3
                d_r = - b * sphi * o3 * ox
                d_b = r * sphi * o3 * ox
                
                Fl_phi = (r2 - br * cphi) * o3 * ox
                pr = 0.0
                pb = 0.0
                
                o = 1.0
                ellippi = cel((sqomm), 1.0 - n, (o), (o))
                ellipe = cel((sqomm), (o), (o), (1.0 - m))
                ellipf = cel((sqomm), (o), (o), (o))
                
                eplusf = alpha * ellipe + beta * ellipf
                eplusf_r = ur * ellipe + vr * ellipf
                eplusf_b = ub * ellipe + vb * ellipf
            else                
                z = sqrt((1.0 - b) * (1.0 + b) - r2 + 2 * br * cphi)
                oz = 1.0 / z
                d = (o3 * (phi * 0.5 - atan(bpr * tphihalf * obmr))
                    - (2 * br * o9) * sphi * z)
                d_r = o9 * b * (-3.0 * ox + 2 * (r2 - br * cphi) * oz - 2 * z) * sphi
                d_b = o9 * r * (3.0 * ox + 2 * (b2 - br * cphi) * oz - 2 * z) * sphi

                Fl_phi = o3 * (1.0 - z ^3) * (r2 - br * cphi) * ox
                    
                pr = (b * sphi * (3.0 - 4 * b2 + b2 * b2 - r2 * (4.0 + r2) + 2 * br * (4.0 - b2 + r2) * cphi)
                   / (9 * y * sqrt(2 * cphi - (b2 + r2 - 1.0) * oy * oy) * (b2 + r2 - 2 * br * cphi)))
                pb = (r * sphi * (-3.0 + 2 * b2 - b2 * b2 + 2 * r2 + r2 * r2 + 2 * br * (-2.0 + b2 - r2) * cphi)
                   / (9 * y * sqrt(2 * cphi - (b2 + r2 - 1.0) * oy * oy) * (b2 + r2 - 2 * br * cphi)))
                                        
                tans = 1.0 / sqrt(m / (sphihalf * sphihalf) - 1.0)
                
                o = 1.0
                ellippi = el3((tans), (sqomm), 1.0 - n)
                ellipe = el2((tans), (sqomm), (o), (1.0 - m))
                ellipf = el2((tans), (sqomm), (o), (o))
                
                eplusf = alpha * ellipe + beta * ellipf
                eplusf_r = ur * ellipe + vr * ellipf
                eplusf_b = ub * ellipe + vb * ellipf

            end
            F_[2] = eplusf + gamma * ellippi + d
            Fl_r = eplusf_r + pr + d_r
            Fl_b = eplusf_b + pb + d_b
        else
        
            z = sqrt((1.0 - bmr) * (1.0 + bmr))
            oz = 1.0 / z
            y = sqrt((1.0 - b) * (1.0 + b) - r2 + 2 * br * cphi)
            oy = 1.0 / y
            x = (b2 + r2 - 2 * br * cphi)
            ox = 1.0 / x
            
            alpha = (7 * r2 + b2 - 4.0) * z * o9
            beta = (r2 * r2 + b2 * b2 + r2 - b2 * (5.0 + 2 * r2) + 1.0) * o9 * oz
            gamma = bpr * o3 * oz * obmr
                        
            n = -4 * br * obmr * obmr
            m = 4 * br * oz * oz
            
            sqomm = sqrt(1.0 - m)
            
            ur = 2 * r * z
            pr = (b * sphi * (3.0 - 4 * b2 + b2 * b2 - r2 * (4.0 + r2) + 2 * br * (4.0 - b2 + r2) * cphi)
               / (9 * y * (1.0 - y) * (1.0 + y)))
                
            ub = z * ((b + 1.0) * (b - 1.0) + r2) / (3 * b)
            vb = (b2 * b2 + ((r - 1.0) * (r + 1.0))^2 - 2 * b2 * (1.0 + r2)) / (3 * b * z)
            pb = r * sphi * (-3.0 + 2 * b2 - b2 * b2 + 2 * r2 + r2 * r2
               + 2 * br * (-2.0 + b2 - r2) * cphi) * ox * o9 * oy 
            
            d = (o3 * (phi * 0.5 - atan(bpr * tphihalf * obmr))
                - 2 * br * o9 * sphi * y)
            d_r = o9 * b * sphi * (-3.0 * ox + 2 * (r2 - br * cphi) * oy - 2 * y)
            d_b = o9 * r * (3.0 * ox + 2 * (b2 - br * cphi) * oy - 2 * y) * sphi
             
            o = 1.0
            ellipe = el2((tphihalf), (sqomm), (o), (1.0 - m))
            ellipf = el2((tphihalf), (sqomm), (o), (o))
            
            if b == r
                Fl_phi = o3 * (1.0 - y^3) * (r2 - br * cphi) * ox
                ellippi = 0.0
                gamma = 0.0
            else
                Fl_phi = o3 * (1.0 - y ^3) * (r2 - br * cphi) * ox
                ellippi = el3((tphihalf), (sqomm), (1.0 - n))
            end
                
            eplusf = alpha * ellipe + beta * ellipf
            eplusf_r = ur * ellipe
            eplusf_b = ub * ellipe + vb * ellipf 
            
            F_[2] = eplusf + gamma * ellippi + d
            Fl_r = eplusf_r + pr + d_r
            Fl_b = eplusf_b + pb + d_b
        end
        
    end
        
    F_rp[2] = Fl_phi * phi_rp
    F_rm[2] = Fl_phi * phi_rm
    F_bp[2] = Fl_phi * phi_bp
    F_bpm[2] = Fl_phi * phi_bpm
    F_theta[2] = Fl_phi * phi_theta
    
    if pflag
        F_bp[3] += Fq_b
        F_rp[3] += Fq_r
        
        F_bp[1] += Fc_b
        F_rp[1] += Fc_r
        
        F_bp[2] += Fl_b
        F_rp[2] += Fl_r
    else
        F_theta[3] += Fq_b * dbm[3]
        F_bpm[3]   += Fq_b * dbm[2]
        F_bp[3]    += Fq_b * dbm[1]
        F_rm[3]    += Fq_r
        
        F_theta[1] += Fc_b * dbm[3]
        F_bpm[1]   += Fc_b * dbm[2]
        F_bp[1]    += Fc_b * dbm[1]
        F_rm[1]    += Fc_r
        
        F_theta[2] += Fl_b * dbm[3]
        F_bpm[2]   += Fl_b * dbm[2]
        F_bp[2]    += Fl_b * dbm[1]
        F_rm[2]    += Fl_r
    end
    
    F[1] = sum(ld .* F_)
    F[2] = sum(ld .* F_rp)
    F[3] = sum(ld .* F_rm)
    F[4] = sum(ld .* F_bp)
    F[5] = sum(ld .* F_bpm)
    F[6] = sum(ld .* F_theta)
    F[7] = - F_[1] + F_[2]
    F[8] = -2 * F_[1] + 2 * F_[2] + F_[3]

    return F
end
