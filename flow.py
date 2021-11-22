if (bp > rp + 1) & (bm > rm + 1):
            ! neither planet nor moon overlap star
            lc(:, i) = f0 * of0
        else if (bpm(i) .gt. rp + rm) then
            if (bp(i) .gt. rp + 1.d0) then
                if (bm(i) .gt. rm + 1.d0) then
                    ! neither planet nor moon overlap star 
                    lc(:, i) = f0 * of0
                else
                    call bm_x(bp(i), bm(i), bpm(i), cth(i), sth(i), dbm)
                    if (bm(i) + rm .le. 1.d0) then
                        ! moon completely overlaps star, planet is outside of star
                        lc(:, i) = (f0 - 2 * Fcomplete(ld, rm, bm(i), dbm, .FALSE.)) * of0
                    else
                        ! moon partially overlaps star, planet is outside of star
                        call kappas_m(rm, bp(i), bm(i), bpm(i), cth(i), sth(i), km, kms, &
                                      km_rm, km_bp, km_bpm, km_theta, &
                                      kms_rm, kms_bp, kms_bpm, kms_theta)
                        lc(:, i) = 2 * (Fstar(ld, pi - kms, 0.d0, -kms_rm, -kms_bp, -kms_bpm, -kms_theta) &
                                        - F(ld, km, rm, bm(i), 0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                            dbm, .FALSE., .TRUE.)) * of0
                    end if
                end if
            else
                if (bm(i) .gt. rm + 1.d0) then
                    if (bp(i) + rp .le. 1.d0) then
                        ! planet completely overlaps star, moon is outside of star
                        lc(:, i) = (f0 - 2 * Fcomplete(ld, rp, bp(i), dbm0, .TRUE.)) * of0
                    else
                        ! planet partially overlaps star, moon is outside of star
                        call kappas_p(rp, bp(i), kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                        lc(:, i) = 2 * (Fstar(ld, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) &
                                        - F(ld, kp, rp, bp(i), kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                            dbm0, .TRUE., .TRUE.)) * of0

                    end if
                else
                    call bm_x(bp(i), bm(i), bpm(i), cth(i), sth(i), dbm)
                    if (bp(i) + rp .le. 1.d0) then
                        if (bm(i) + rm .le. 1.d0) then
                            ! moon and planet both completely overlap star, they do not overlap each othe
                            lc(:, i) = (f0 - 2 * (Fcomplete(ld, rm, bm(i), dbm, .FALSE.) &
                                  + Fcomplete(ld, rp, bp(i), dbm0, .TRUE.))) * of0
                        else
                            ! planet completely overlaps star, moon partially overlaps star, they do not overlap each other
                            call kappas_m(rm, bp(i), bm(i), bpm(i), cth(i), sth(i), km, kms, &
                                      km_rm, km_bp, km_bpm, km_theta, &
                                      kms_rm, kms_bp, kms_bpm, kms_theta)
                            lc(:, i) = 2 * (Fstar(ld, pi - kms, 0.d0, -kms_rm, -kms_bp, -kms_bpm, -kms_theta) &
                                            - F(ld, km, rm, bm(i), 0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                                dbm, .FALSE., .TRUE.) &
                                            - Fcomplete(ld, rp, bp(i), dbm0, .TRUE.)) * of0
                        end if
                    else
                        if (bm(i) + rm .le. 1.d0) then
                            ! planet partially overlaps star, moon fully overlaps star, they do not overlap each other
                            call kappas_p(rp, bp(i), kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                            lc(:, i) = 2 * (Fstar(ld, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) & 
                                            - F(ld, kp, rp, bp(i), kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                dbm0, .TRUE., .TRUE.) &
                                            - Fcomplete(ld, rm, bm(i), dbm, .FALSE.)) * of0
                        else
                            ! moon and planet both partially overlap star, but not each other
                            call kappas_p(rp, bp(i), kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                            call kappas_m(rm, bp(i), bm(i), bpm(i), cth(i), sth(i), km, kms, &
                                      km_rm, km_bp, km_bpm, km_theta, &
                                      kms_rm, kms_bp, kms_bpm, kms_theta)
                            lc(:, i) = 2 * (Fstar(ld, pi - (kps + kms), -kps_rp, -kms_rm, &
                                                  -(kps_bp + kms_bp), -kms_bpm, -kms_theta) &
                                            - F(ld, kp, rp, bp(i), kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                dbm0, .TRUE., .TRUE.) &
                                            - F(ld, km, rm, bm(i), 0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                                dbm, .FALSE., .TRUE.)) * of0
                        end if
                    end if
                end if
            end if
        else
            if (bp(i) .gt. rp + 1.d0) then
                if (bm(i) .gt. rm + 1.d0) then
                    ! neither moon nor planet overlap star
                    lc(:, i) = f0
                else
                    ! moon partially overlaps star, planet does not overlap star
                    call bm_x(bp(i), bm(i), bpm(i), cth(i), sth(i), dbm)
                    call kappas_m(rm, bp(i), bm(i), bpm(i), cth(i), sth(i), km, kms, &
                                      km_rm, km_bp, km_bpm, km_theta, &
                                      kms_rm, kms_bp, kms_bpm, kms_theta)
                    lc(:, i) = 2 * (Fstar(ld, pi - kms, 0.d0, -kms_rm, -kms_bp, -kms_bpm, -kms_theta) &
                                    - F(ld, km, rm, bm(i), 0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                            dbm, .FALSE., .TRUE.)) * of0

                end if
            else
                if (bm(i) .gt. rm + 1.d0) then
                    if (bp(i) + rp .le. 1.d0) then
                        ! planet fully overlaps star, moon does not overlap star
                        lc(:, i) = (f0 - 2 * Fcomplete(ld, rp, bp(i), dbm0, .TRUE.)) * of0
                    else
                        ! planet partially overlaps star, moon is outside of star
                        call kappas_p(rp, bp(i), kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                        lc(:, i) = 2 * (Fstar(ld, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) &
                                        - F(ld, kp, rp, bp(i), kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                            dbm0, .TRUE., .TRUE.)) * of0
                    end if
                else
                    if (bp(i) + rp .le. 1.d0) then
                        if (bm(i) + rm .le. 1.d0) then
                            if (bpm(i) + rm .le. rp) then
                                ! moon and planet both overlap star, moon fully overlapped by planet
                                lc(:, i) = (f0 - 2 * Fcomplete(ld, rp, bp(i), dbm0, .TRUE.)) * of0
                            else
                                ! moon and planet both overlap star, moon and planet partially overlap each other 
                                ! bookmark
                                call bm_x(bp(i), bm(i), bpm(i), cth(i), sth(i), dbm)
                                call phis(rp, rm, bp(i), bm(i), bpm(i), cth(i), sth(i), pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, &
                                          pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta)
                                lc(:, i) = (f0 - Arc(ld, pp1, pp2, rp, bp(i), &
                                                     pp_rp, pp_rm, 0.d0, pp_bpm, 1.d0, &
                                                     -pp_rp, -pp_rm, 0.d0, -pp_bpm, 1.d0, &
                                                     dbm0, .TRUE., .FALSE., .FALSE.) &
                                                - Arc(ld, pm1, pm2, rm, bm(i), &
                                                     pm_rp, pm_rm, thetam_bp, pm_bpm + thetam_bpm, thetam_theta, &
                                                     -pm_rp, -pm_rm, thetam_bp, -pm_bpm + thetam_bpm, thetam_theta, &
                                                     dbm, .FALSE., .FALSE., .FALSE.)) * of0
                            end if
                        else
                            call bm_x(bp(i), bm(i), bpm(i), cth(i), sth(i), dbm)
                            call phis(rp, rm, bp(i), bm(i), bpm(i), cth(i), sth(i), pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, &
                                      pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta)
                            call kappas_m(rm, bp(i), bm(i), bpm(i), cth(i), sth(i), km, kms, &
                                      km_rm, km_bp, km_bpm, km_theta, &
                                      kms_rm, kms_bp, kms_bpm, kms_theta)
                            lc(:, i) = (2 * Fstar(ld, pi - kms, 0.d0, -kms_rm, -kms_bp, -kms_bpm, -kms_theta) &
                                        - Arc(ld, -km, pm2, rm, bm(i), &
                                              0.d0, -km_rm, -km_bp, -km_bpm, -km_theta, &
                                              -pm_rp, -pm_rm, thetam_bp, -pm_bpm + thetam_bpm, thetam_theta, &
                                              dbm, .FALSE., .TRUE., .FALSE.) &
                                        - Arc(ld, pm1, km, rm, bm(i), &
                                              pm_rp, pm_rm, thetam_bp, pm_bpm + thetam_bpm, thetam_theta, &
                                              0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                              dbm, .FALSE., .FALSE., .TRUE.) &
                                        - Arc(ld, pp1, pp2, rp, bp(i), &
                                              pp_rp, pp_rm, 0.d0, pp_bpm, 1.d0, &
                                              -pp_rp, -pp_rm, 0.d0, -pp_bpm, 1.d0, &
                                              dbm0, .TRUE., .FALSE., .FALSE.)) * of0
                        end if
                    else
                        if (bm(i) + rm .le. 1.d0) then 
                            if (bpm(i) + rm .le. rp) then
                                ! planet partially overlaps star, moon fully overlaps star but is completely overlapped by planet 
                                call kappas_p(rp, bp(i), kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                                lc(:, i) = 2 * (Fstar(ld, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) &
                                                - F(ld, kp, rp, bp(i), kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                    dbm0, .TRUE., .TRUE.)) * of0
                            else
                                ! planet partially overlaps star, moon fully overlaps star and only partially overlaps planet
                                call bm_x(bp(i), bm(i), bpm(i), cth(i), sth(i), dbm)
                                call phis(rp, rm, bp(i), bm(i), bpm(i), cth(i), sth(i), pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, &
                                      pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta)
                                call kappas_p(rp, bp(i), kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                                lc(:, i) = (2 * Fstar(ld, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) &
                                            - Arc(ld, -kp, pp2, rp, bp(i), &
                                                  -kp_rp, 0.d0, -kp_bp, 0.d0, 0.d0, &
                                                  -pp_rp, -pp_rm, 0.d0, -pp_bpm, 1.d0, &
                                                  dbm0, .TRUE., .TRUE., .FALSE.) &
                                            - Arc(ld, pp1, kp, rp, bp(i), &
                                                  pp_rp, pp_rm, 0.d0, pp_bpm, 1.d0, &
                                                  kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                  dbm0, .TRUE., .FALSE., .TRUE.) &
                                            - Arc(ld, pm1, pm2, rm, bm(i), &
                                                  pm_rp, pm_rm, thetam_bp, pm_bpm + thetam_bpm, thetam_theta, &
                                                  -pm_rp, -pm_rm, thetam_bp, -pm_bpm + thetam_bpm, thetam_theta, &
                                                  dbm, .FALSE., .FALSE., .FALSE.)) * of0
                            end if
                        else
                            if (bpm(i) + rm .le. rp) then
                                ! planet and moon both partially overlap star but moon is fully overlapped by the planet
                                call kappas_p(rp, bp(i), kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                                lc(:, i) = 2 * (Fstar(ld, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) &
                                                - F(ld, kp, rp, bp(i), kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                    dbm0, .TRUE., .TRUE.)) * of0
                            else
                                !call compute_theta(rp,  bp(i), theta, phip, theta_bp, theta_rp, phip_bp, phip_rp)
                                !call compute_theta(rm,  bm(i), theta, phim, theta_bm, theta_rm, phim_bm, phim_rm)
                                call kappas_p(rp, bp(i), kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                                call kappas_m(rm, bp(i), bm(i), bpm(i), cth(i), sth(i), km, kms, &
                                      km_rm, km_bp, km_bpm, km_theta, &
                                      kms_rm, kms_bp, kms_bpm, kms_theta)
                                
                                a = bm(i)
                                b = bp(i)
                                c = bpm(i)
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
                                phi = Atan2(delta, (bm(i) - bpm(i)) * (bm(i) + bpm(i)) + bp(i) * bp(i))                                
                                
                                call phis(rp, rm, bp(i), bm(i), bpm(i), cth(i), sth(i), pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, &
                                          pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta)
                                
                                if (phi + kms .le. kps) then  
                                        if (pp2 .gt. kp) then
                                            ! planet and moon both partially overlap the star and each other but the 
                                            ! moon-star overlap is contained within the planet-star overlap
                                            lc(:, i) = 2 * (Fstar(ld, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) &
                                                            - F(ld, kp, rp, bp(i), kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                                dbm0, .TRUE., .TRUE.)) * of0
                                        else
                                            ! planet and moon both partially overlap star and each other but the 
                                            ! planet-star intersections are overlapped by the planet
                                            call bm_x(bp(i), bm(i), bpm(i), cth(i), sth(i), dbm)
                                            lc(:, i) = (2 * Fstar(ld, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) &
                                                        - Arc(ld, -kp, pp2, rp, bp(i), &
                                                              -kp_rp, 0.d0, -kp_bp, 0.d0, 0.d0, &
                                                              -pp_rp, -pp_rm, 0.d0, -pp_bpm, 1.d0, &
                                                              dbm0, .TRUE., .TRUE., .FALSE.) &
                                                        - Arc(ld, pp1, kp, rp, bp(i), &
                                                              pp_rp, pp_rm, 0.d0, pp_bpm, 1.d0, &
                                                              kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                              dbm0, .TRUE., .FALSE., .TRUE.) &
                                                        - Arc(ld, pm1, pm2, rm, bm(i), &
                                                              pm_rp, pm_rm, thetam_bp, pm_bpm + thetam_bpm, thetam_theta, &
                                                              -pm_rp, -pm_rm, thetam_bp, -pm_bpm + thetam_bpm, thetam_theta, &
                                                              dbm, .FALSE., .FALSE., .FALSE.)) * of0
                                        end if
                                else if (phi + kps .le. kms) then
                                    call bm_x(bp(i), bm(i), bpm(i), cth(i), sth(i), dbm)
                                    if ((bp(i) - rp) .le. (bm(i) - rm)) then
                                        ! planet and moon both partially overlap the star and each other but the 
                                        ! planet-star intersections are overlapped by the moon
                                        lc(:, i) = (2 * Fstar(ld, pi - kms, 0.d0, -kms_rm, -kms_bp, -kms_bpm, -kms_theta) &
                                                        - Arc(ld, -km, pm2, rm, bm(i), &
                                                              0.d0, -km_rm, -km_bp, -km_bpm, -km_theta, &
                                                              -pm_rp, -pm_rm, thetam_bp, -pm_bpm + thetam_bpm, thetam_theta, &
                                                              dbm, .FALSE., .TRUE., .FALSE.) &
                                                        - Arc(ld, pm1, km, rm, bm(i), &
                                                              pm_rp, pm_rm, thetam_bp, pm_bpm + thetam_bpm, thetam_theta, &
                                                              0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                                              dbm, .FALSE., .FALSE., .TRUE.) &
                                                        - Arc(ld, pp1, pp2, rp, bp(i), &
                                                              pp_rp, pp_rm, 0.d0, pp_bpm, 1.d0, &
                                                              -pp_rp, -pp_rm, 0.d0, -pp_bpm, 1.d0, &
                                                              dbm0, .TRUE., .FALSE., .FALSE.)) * of0
                                    else
                                        ! planet and moon both partially overlap the star and each other but 
                                        ! the planet-star overlap is  entirely within the moon-star overlap
                                        lc(:, i) = 2 * (Fstar(ld, pi - kms, 0.d0, -kms_rm, -kms_bp, -kms_bpm, -kms_theta) &
                                                        - F(ld, km, rm, bm(i), 0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                                            dbm, .FALSE., .TRUE.)) * of0
                                    end if
                                else
                                    ! bookmark
                                    d1 = rm * rm + bm(i) * bm(i) - 2 * rm * bm(i) * Cos(pm2)
                                    d2 = rm * rm + bm(i) * bm(i) - 2 * rm * bm(i) * Cos(pm1)
                                    call bm_x(bp(i), bm(i), bpm(i), cth(i), sth(i), dbm)
                                    if ((d1 .gt. 1.d0) .AND. (d2 .gt. 1.d0)) then
                                        ! planet and moon both partially overlap star and each other, 
                                        ! but the planet/moon overlap does not overlap the star
                                        lc(:, i) = 2 * (Fstar(ld, pi - (kps + kms), -kps_rp, -kms_rm, &
                                                              -(kps_bp + kms_bp), -kms_bpm, -kms_theta) &
                                                        - F(ld, kp, rp, bp(i), kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                            dbm0, .TRUE., .TRUE.) &
                                                        - F(ld, km, rm, bm(i), 0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                                            dbm0, .FALSE., .TRUE.)) * of0
                                    else if ((d1 .le. 1.d0) .AND. (d2 .le. 1.d0)) then
                                        ! planet and moon both partially overlap star and each other, 
                                        ! with the planet/moon overlap fully overlapping the star
                                        lc(:, i) = (2 * Fstar(ld, pi - (kps + kms), -kps_rp, -kms_rm, &
                                                              -(kps_bp + kms_bp), -kms_bpm, -kms_theta) &
                                                        - Arc(ld, -km, -pm1, rm, bm(i), &
                                                              0.d0, -km_rm, -km_bp, -km_bpm, -km_theta, &
                                                              -pm_rp, -pm_rm, -thetam_bp, -pm_bpm - thetam_bpm, -thetam_theta, &
                                                              dbm, .FALSE., .TRUE., .FALSE.) &
                                                        - Arc(ld, -pm2, km, rm, bm(i), &
                                                              pm_rp, pm_rm, -thetam_bp, pm_bpm - thetam_bpm, -thetam_theta, &
                                                              0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                                              dbm, .FALSE., .FALSE., .TRUE.) &
                                                        - Arc(ld, pp1, kp, rp, bp(i), &
                                                              pp_rp, pp_rm, 0.d0, pp_bpm, 1.d0, &
                                                              kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                              dbm0, .TRUE., .FALSE., .TRUE.) &
                                                        - Arc(ld, -kp, pp2, rp, bp(i), &
                                                              -kp_rp, 0.d0, -kp_bp, 0.d0, 0.d0, &
                                                              -pp_rp, -pp_rm, 0.d0, -pp_bpm, 1.d0, &
                                                              dbm0, .TRUE., .TRUE., .FALSE.)) * of0
                                    else
                                        ! planet and moon both partially overlap star and each other, 
                                        ! with the planet/moon overlap partially overlapping the star
                                        
                                        ! there might be a mistake somewhere in here... 
                                        phi_bpm = bp(i) * sth(i) / (bm(i) * bm(i))
                                        phi_theta = bpm(i) * (bp(i) * cth(i) - bpm(i)) / (bm(i) * bm(i))
                                        phi_bp = - bpm(i) * sth(i) / (bm(i) * bm(i))
                                        
                                        lc(:, i) = (2 * Fstar(ld, pi - 0.5 * (kps + kms + phi), &
                                                              -0.5 * kps_rp, -0.5 * kms_rm, &
                                                              -0.5 * (kps_bp + kms_bp + phi_bp), -0.5 * (kms_bpm + phi_bpm), &
                                                              -0.5 * (kms_theta + phi_theta)) &
                                                        - Arc(ld, -pm2, km, rm, bm(i), &
                                                              pm_rp, pm_rm, -thetam_bp, pm_bpm - thetam_bpm, -thetam_theta, &
                                                              0.d0, km_rm, km_bp, km_bpm, km_theta,  &
                                                              dbm, .FALSE., .FALSE., .TRUE.) &
                                                        - Arc(ld, -kp, pp2, rp, bp(i), &
                                                              -kp_rp, 0.d0, -kp_bp, 0.d0, 0.d0, &
                                                              -pp_rp, -pp_rm, 0.d0, -pp_bpm, 1.d0, &
                                                              dbm0, .TRUE., .TRUE., .FALSE.)) * of0
                                    end if
                                end if
                            end if
                        end if
                    end if
                end if
            end if  
        end if