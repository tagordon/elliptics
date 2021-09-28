module kepler
use iso_c_binding
implicit none

real*8, parameter :: pi = 3.14159265358979323846, G = 8.8876413d-10
real*8, parameter :: o3 = 0.33333333333333333333
real*8, parameter :: days_in_year = 365.256, earths_in_sun = 332946.08

contains

subroutine solve_kepler(t, n, t0, ecc, a, r, cosf, sinf, j)

    integer :: j, i
    real*8, dimension(j) :: t, M, E, sE, cE
    real*8, dimension(j), intent(out) :: r, cosf, sinf
    real*8, dimension(j) :: y, y2, denom
    real*8 :: n, t0, ecc, a
    real*8 :: tol, x, err
    
    M = n * (t - t0)
    if (ecc .lt. 1.d-5) then
        cosf = Cos(M)
        sinf = Sin(M)
        r = a
    else
        tol = 1.d-7
        
        E = M + ecc
        x = Sqrt((1.d0 + ecc) / (1.d0 - ecc))
        
        do i=1,j,1
            err = 1.d0
            do while (err .gt. tol)
                err = - (E(i) - ecc * Sin(E(i)) - M(i)) / (1.d0 - ecc * Cos(E(i)))
                E(i) = E(i) + err
            end do
        end do    
        
        cE = Cos(E)
        y = x * Sin(E) / (1.d0 + cE)
        y2 = y * y
        denom = 1.d0 / (y2 + 1.d0)
        cosf = (1.d0 - y2) * denom
        sinf = 2 * y * denom
        r = a * (1.d0 - ecc * cE)
    end if
end

subroutine solve_kepler_markley(t, n, t0, ecc, a, rr, sinf, cosf, j)

    integer :: j, i
    real*8 :: t(j), rr(j), M(j), E(j), err(j), sinf(j), cosf(j)
    real*8 :: M2(j), alpha(j), d(j), alphad(j), r(j), q(j), q2(j), w(j), dE(j)
    real*8 :: f0(j), f1(j), f2(j), f3(j), d3(j), d4(j), d42(j), sE(j), cE(j), y(j), y2(j), denom(j)
    real*8 :: n, t0, ecc, a
    real*8 :: fact1, fact2, ome, x

    M = n * (t - t0)
    if (ecc .lt. 1.d-5) then
        cosf = Cos(M)
        sinf = Sin(M)
        rr = a
    else
        ome = 1.d0 - ecc
        fact1 = 3 * pi / (pi - 6 / pi)
        fact2 = 1.6 / (pi - 6 / pi)
        x = Sqrt((1.d0 + ecc) / (1.d0 - ecc))
        
        M2 = M * M
        alpha = fact1 + fact2 * (pi - M) / (1.d0 + ecc)
        d = 3 * ome + alpha * ecc
        alphad = alpha * d
        r = (3 * alphad * (d - ome) + M2) * M
        q = 2 * alphad * ome - M2
        q2 = q * q
        w = (Abs(r) + Sqrt(q2 * q + r * r)) ** (2.d0 / 3.d0)
        E = (2 * r * w / (w * w + w * q + q2) + M) / d
        
        sE = E - Sin(E)
        cE = 1.d0 - Cos(E)
        f0 = ecc * sE + E * ome - M
        f1 = ecc * cE + ome
        f2 = ecc * (E - sE)
        f3 = 1.d0 - f1
        d3 = -f0 / (f1 - 0.5 * f0 * f2 / f1)
        d4 = -f0 / (f1 + 0.5 * d3 * f2 + (d3 * d3) * f3 / 6.d0)
        d42 = d4 * d4
        dE = -f0 / (f1 + 0.5 * d4 * f2 + d4 * d4 * f3 / 6.d0 - d42 * d4 * f2 / 24.d0)
        E = E + dE
            
        sE = Sin(E)
        cE = Cos(E)
        y = x * sE / (1 + cE)
        y2 = y * y
        denom = 1.d0 / (y2 + 1.d0)
        cosf = (1.d0 - y2) * denom
        sinf = 2 * y * denom
        rr = a * (1.d0 - ecc * cE)
    end if

end

real*8 function tt(ecc, w, p, t0, n)

    real*8 :: ecc, w, p, t0, n
    real*8 :: x, E, M
    
    if (ecc .lt. 1.d-5) then
        tt = p * (pi * 0.5 - w) / (2 * pi) + t0
        return
    else
        x = (1.d0 - (1.d0 - ecc**2.d0)) / (ecc * (1.d0 + ecc * Cos(pi * 0.5 - w)))
        if (x .gt. 1.d0) then
            x = 1.d0
        end if
        E = Acos(x)
        M = E - ecc * Sin(E)
        tt = M / n + t0
        return
    end if
end 

real*8 function dt(w, ep, em, rad, ap, am, ms, mp, mm)

    real*8 :: w, ep, em, rad, ap, am, ms, mp, mm
    real*8 :: ft, rt, df, fa, dfa, dfb, h
    
    ft = 0.5 * pi - w
    rt = ap * (1.d0 - ep**2.d0) / (1.d0 + ep * Cos(ft))
    df = Asin((rad + (em + 1.d0) * am) / rt)
    fa = ft - df
    dfa = 2 * Sqrt(1.d0 - ep**2.d0) * Atan(Sqrt((1.d0 - ep) / (1.d0 + ep)) * Tan(fa * 0.5)) - &
        (ep * (1.d0 - ep ** 2.d0) * Sin(fa)) / (1.d0 + ep * Cos(fa))
    dfb = 2 * Sqrt(1.d0 - ep**2.d0) * Atan(sqrt((1.d0 - ep) / (1.d0 + ep)) * Tan(ft * 0.5)) - &
        (ep * (1.d0 - ep**2.d0) * Sin(ft)) / (1.d0 + ep * Cos(ft))
    h = Sqrt(G * (ms + mp + mm) * ap * (1.d0 - ep**2.d0))
    dt = (dfb - dfa) * ap**2.d0 / h

end

subroutine pscoords(t, n, t0, e, a, mprim, msec, w, omega, i, xp, yp, zp, xs, ys, zs, j)

    integer :: j
    real*8 :: t(j), xp(j), yp(j), zp(j), xs(j), ys(j), zs(j), x(j), y(j), z(j)
    real*8 :: n, t0, e, a, mprim, msec, w, omega, i, mrprim, mrsec, comega, somega, ci, cw, sw
    real*8 :: cosf(j), sinf(j), r(j), f(j), cosfw(j), sinfw(j)
    
    call solve_kepler(t, n, t0, e, a, r, cosf, sinf, j)
    
    mrprim = -msec / (mprim + msec)
    mrsec = mprim / (mprim + msec)
    
    comega = Cos(omega)
    somega = Sin(omega)
    ci = Cos(i)
    cw = Cos(w)
    sw = Sin(w)
    
    cosfw = (cw * cosf - sw * sinf)
    sinfw = (sw * cosf + sinf * cw)
    
    x = -r * (comega * cosfw - somega * sinfw * ci)
    y = -r * (somega * cosfw + comega * sinfw * ci)
    z = r * sinfw * Sin(i)
    
    xp = mrprim * x
    yp = mrprim * y
    zp = mrprim * z
    
    xs = mrsec * x
    ys = mrsec * y
    zs = mrsec * z

end

subroutine coords_no_z(t, n, t0, e, a, mprim, msec, w, omega, i, xp, yp, xs, ys, j)

    integer :: j
    real*8 :: t(j), xp(j), yp(j), zp(j), xs(j), ys(j), zs(j), x(j), y(j), z(j)
    real*8 :: n, t0, e, a, mprim, msec, w, omega, i, mrprim, mrsec, comega, somega, ci, cw, sw
    real*8 :: cosf(j), sinf(j), r(j), f(j), cosfw(j), sinfw(j)
    
    call solve_kepler(t, n, t0, e, a, r, cosf, sinf, j)
    
    mrprim = -msec / (mprim + msec)
    mrsec = mprim / (mprim + msec)
    
    comega = Cos(omega)
    somega = Sin(omega)
    ci = Cos(i)
    cw = Cos(w)
    sw = Sin(w)
    
    cosfw = (cw * cosf - sw * sinf)
    sinfw = (sw * cosf + sinf * cw)
    
    x = -r * (comega * cosfw - somega * sinfw * ci)
    y = -r * (somega * cosfw + comega * sinfw * ci)
    
    xp = mrprim * x
    yp = mrprim * y
    
    xs = mrsec * x
    ys = mrsec * y
    
end

subroutine pscoords_no_z(t, n, t0, e, a, mprim, msec, w, omega, i, xp, yp, xs, ys, j)

    integer :: j
    real*8 :: t(j), xp(j), yp(j), xs(j), ys(j), x(j), y(j)
    real*8 :: n, t0, e, a, mprim, msec, w, omega, i, mrprim, mrsec, comega, somega, cw, sw, ci
    real*8 :: cosf(j), sinf(j), r(j), f(j), cosfw(j), sinfw(j)
    
    call solve_kepler(t, n, t0, e, a, r, cosf, sinf, j)
    
    mrprim = -msec / (mprim + msec)
    mrsec = mprim / (mprim + msec)
    
    comega = Cos(omega)
    somega = Sin(omega)
    cw = Cos(w)
    sw = Sin(w)
    ci = Cos(i)
    
    cosfw = (cw * cosf - sw * sinf)
    sinfw = (sw * cosf + sinf * cw)
    
    x = -r * (comega * cosfw - somega * sinfw * ci)
    y = -r * (somega * cosfw + comega * sinfw * ci)
    
    xp = mrprim * x
    yp = mrprim * y
    
    xs = mrsec * x
    ys = mrsec * y
end

subroutine spcoords_no_z(t, n, t0, e, a, mprim, msec, w, omega, i, j)

    integer :: j
    real*8 :: t(j), xs(j), ys(j), x(j), y(j)
    real*8 :: n, t0, e, a, mprim, msec, w, omega, i, mrprim, mrsec, comega, somega, cw, sw, ci
    real*8 :: cosf(j), sinf(j), r(j), f(j), cosfw(j), sinfw(j)
    
    call solve_kepler(t, n, t0, e, a, r, cosf, sinf, j)
    
    mrprim = -msec / (mprim + msec)
    mrsec = mprim / (mprim + msec)
    
    comega = Cos(omega)
    somega = Sin(omega)
    cw = Cos(w)
    sw = Sin(w)
    ci = Cos(i)
    
    cosfw = (cw * cosf - sw * sinf)
    sinfw = (sw * cosf + sinf * cw)
    
    x = -r * (comega * cosfw - somega * sinfw * ci)
    y = -r * (somega * cosfw + comega * sinfw * ci)
    
    xs = mrsec * x
    ys = mrsec * y
end

subroutine get_coords(t, ms, t0p, ep, Pp, Op, wp, ip, mp, &
    &t0m, em, Pm, wm, Om, im, mm, j, &
    &xs, ys, zs, xp, yp, zp, xm, ym, zm) bind(C, name="coords")

    integer (c_int), bind(C) :: j
    real (c_double), bind(C) :: ms, t0p, ep, Pp, Op, wp, ip, mp 
    real (c_double), bind(C) :: t0m, em, Pm, Om, wm, im, mm
    real (c_double), bind(C) :: t(j)
    real*8 :: np, nm, ap, am
    real*8 :: xbc(j), ybc(j), zbc(j)
    real (c_double), bind(C), intent(out) :: xs(j), ys(j), zs(j), xp(j), yp(j), zp(j), xm(j), ym(j), zm(j)
    
    np = 2 * pi / Pp
    nm = 2 * pi / Pm
    ap = (G * (ms + mp) / (np ** 2)) ** o3
    am = (G * (mp + mm) / (nm ** 2)) ** o3
        
    call pscoords(t, np, t0p, ep, ap, ms, mp, wp, Op, ip, xs, ys, zs, xbc, ybc, zbc, j)
    call pscoords(t, nm, t0m, em, am, mp, mm, wm, Om, im, xp, yp, zp, xm, ym, zm, j)
    
    xp = xbc + xp
    yp = ybc + yp
    zp = zbc + zp
    
    xm = xbc + xm
    ym = ybc + ym
    zm = zbc + zm
    
1 end

subroutine get_impacts(t, ms, t0p, ep, Pp, Op, wp, ip, mp, &
    &t0m, em, Pm, wm, Om, im, mm, j, x, y, xbc, ybc, bp2, bm2, bpm2) bind(C, name="impacts")

    integer (c_int), bind(C) :: j
    real (c_double), bind(C) :: ms, t0p, ep, Pp, Op, wp, ip, mp 
    real (c_double), bind(C) :: t0m, em, Pm, Om, wm, im, mm
    real (c_double), bind(C) :: t(j)
    real (c_double), bind(C), intent(out), dimension(j) :: bp2, bm2, bpm2
    real*8 :: np, nm, ap, am
    real (c_double), bind(C), intent(out), dimension(j) :: xbc, ybc, x, y
    real*8, dimension(j) :: cosfbc, sinfbc, rbc, cosfwbc, sinfwbc
    real*8, dimension(j) :: cosfm, sinfm, rm, cosfwm, sinfwm
    real*8 :: mrp, mrm, comegap, somegap, cwp, swp, cip, comegam, somegam, cwm, swm, cim
    
    np = 2 * pi / Pp
    nm = 2 * pi / Pm
    ap = (G * (ms + mp + mm) / (np ** 2)) ** o3
    am = (G * (mp + mm) / (nm ** 2)) ** o3
    
    comegap = Cos(Op)
    somegap = Sin(Op)
    cwp = Cos(wp)
    swp = Sin(wp)
    cip = Cos(ip)
    
    comegam = Cos(Om)
    somegam = Sin(Om)
    cwm = Cos(wm)
    swm = Sin(wm)
    cim = Cos(im)
    
    mrp = -mm / (mp + mm)
    mrm = mp / (mp + mm)
    
    call solve_kepler(t, np, t0p, ep, ap, rbc, cosfbc, sinfbc, j)
    
    cosfwbc = cwp * cosfbc - swp * sinfbc
    sinfwbc = swp * cosfbc + sinfbc * cwp
    xbc = -rbc * (comegap * cosfwbc - somegap * sinfwbc * cip)
    ybc = -rbc * (somegap * cosfwbc + comegap * sinfwbc * cip)
    
    call solve_kepler(t, nm, t0m, em, am, rm, cosfm, sinfm, j)
    
    cosfwm = cwm * cosfm - swm * sinfm
    sinfwm = swm * cosfm + sinfm * cwm
    x = -rm * (comegam * cosfwm - somegam * sinfwm * cim)
    y = -rm * (somegam * cosfwm + comegam * sinfwm * cim)
    
    bp2 = (xbc + mrm * x)**2 + (ybc + mrm * y)**2
    bm2 = (xbc + mrp * x)**2 + (ybc + mrp * y)**2
    bpm2 = x**2 + y**2
    
end

end