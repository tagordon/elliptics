! These routines might be useful for implementing a scipy version of the bulirsch elliptic integrals
! Use ellip.f90 instead sfor the photodynamics code

module bulirsch
use iso_c_binding
use ellip
use phot
implicit none

contains

subroutine test_integral(phi1, phi2, r, b, res, j) bind(C, name="test_int")
    
    integer (c_int), bind(C) :: j
    integer :: i
    real (c_double), bind(C) :: phi1(j), phi2(j), r(j), b(j)
    real (c_double), bind(c), dimension(j), intent(out) :: res
    
    do i=1,j,1
        res(i) = F(phi2(i), r(i), b(i)) + F(phi1(i), r(i), b(i))
    end do
    return
end

subroutine f_bulirsch(phi, m, e, j) bind(C, name="f_burl")

    integer (c_int), bind(C) :: j
    integer :: i
    real (c_double), bind(C) :: phi(j), m(j)
    real (c_double), bind(C), dimension(j), intent(out) :: e
    
    do i=1,j,1
        e(i) = el1(tan(phi(i)), sqrt(1.d0 - m(i)))
    end do
    !call el1(tan(phi), sqrt(1 - m), e, j)
    return
end

subroutine e_bulirsch(phi, m, e, j) bind(C, name="e_burl")

    integer (c_int), bind(C) :: j
    integer :: i
    real (c_double), bind(C) :: phi(j), m(j)
    real*8 :: o(j)
    real (c_double), bind(C), dimension(j), intent(out) :: e
    
    o = 1.d0
    do i=1,j,1
        e(i) = el2(tan(phi(i)), sqrt(1.d0 - m(i)), o(i), 1.d0 - m(i))
    end do
    !call el2(tan(phi), sqrt(1 - m), o, 1 - m, e, j)
    return
end

subroutine p_bulirsch(phi, n, m, e, j) bind(C, name="p_burl")

    integer (c_int), bind(C) :: j
    integer :: i
    real (c_double), bind(C) :: phi(j), m(j), n(j)
    real (c_double), bind(C), dimension(j), intent(out) :: e
    
    do i=1,j,1
        e(i) = el3(tan(phi(i)), sqrt(1.d0 - m(i)), 1.d0 - n(i))
    end do
    !call el3(tan(phi), sqrt(1 - m), 1 - n, e, j)
    return
end

subroutine el1_tmp(x, kc, e, j) bind(C, name="el1")

    integer (c_int), bind(C) :: j
    real (c_double), bind(C) :: x(j), kc(j)
    real (c_double), bind(C), dimension(j), intent(out) :: e
    
    real*8 :: g(j), m(j), y(j)
    real*8 :: yi, ei, mi, kci, gi, xi
    real*8 :: D = 8
    real*8 :: ca, cb
    real*8 :: pi = 3.1415926535897932384626433832795
    integer :: l, i
    
    ca = 10**(-D/2.d0)
    cb = 10**(-D+2.d0)
    
    y = abs(1.d0/x)
    kc = abs(kc)
    m = 1.d0
    
    do i=1,j,1
    
        if (kc(i) == 0.d0) then
            e(i) = log(x(i) + 1.d0 / cos(atan(x(i))))
            goto 2
        end if
    
        yi = y(i)
        ei = e(i)
        mi = m(i)
        kci = kc(i)
        gi = g(i)
        xi = x(i)
        l = 0
        
1       ei = mi * kci
        gi = mi
        mi = kci + mi
        yi = - (ei/yi) + yi
        if (yi == 0) then
            yi = sqrt(ei) * cb
        end if 
        if (abs(gi - kci) .gt. ca * gi) then
            kci = 2 * sqrt(ei)
            l = 2 * l
            if (yi .lt. 0) then
                l = l + 1
            end if 
            goto 1
        end if
        if (yi .lt. 0) then 
            l = l + 1
        end if 
        ei = (atan(mi/yi) + pi * l) / mi
        if (xi .lt. 0) then
            e(i) = -ei
        else
            e(i) = ei
        end if
2   end do
    
    return
end

subroutine el2_tmp(x, kc, a, b, e, j) bind(C, name="el2")

    integer (c_int), bind(C) :: j
    real (c_double), bind(C) :: x(j), kc(j), a(j), b(j)
    real (c_double), bind(C), dimension(j), intent(out) :: e
    
    real*8 :: c(j), dd(j), f(j), g(j), ik(j), m(j), p(j), y(j), z(j)
    real*8 :: ai, bi, ci, di, fi, iki, ppi, yi, zi, ei, mi, kci, gi, xi
    real*8 :: D = 8
    real*8 :: ca, cb
    real*8 :: pi = 3.1415926535897932384626433832795
    integer :: l, i
    
    ca = 10**(-D/2.d0)
    cb = 10**(-D+2.d0)
    
    c = x**2
    dd = c + 1.d0
    p = sqrt((1.d0 + kc**2 * c)/dd)
    dd = x / dd
    c = dd / (2.d0 * p)
    z = a - b
    ik = a
    a = (b + a) / 2.d0
    y = abs(1.d0/x)
    f = 0.d0
    kc = abs(kc)
    m = 1.d0
    
    do i=1,j,1
    
        if (kc(i) == 0.d0) then
            e(i) = sin(atan(x(i)))
            goto 2
        end if
    
        ai = a(i)
        bi = b(i)
        ci = c(i)
        di = dd(i)
        ppi = p(i)
        zi = z(i)
        iki = ik(i)
        fi = f(i)
        yi = y(i)
        ei = e(i)
        mi = m(i)
        kci = kc(i)
        gi = g(i)
        xi = x(i)
        l = 0
        
1       bi = iki * kci + bi
        ei = mi * kci
        gi = ei / ppi
        di = fi * gi + di
        fi = ci
        iki = ai
        ppi = gi + ppi
        ci = (di / ppi + ci) / 2.d0
        gi = mi
        mi = kci + mi
        ai = (bi / mi + ai) / 2.d0
        yi = - (ei/yi) + yi
        if (yi == 0) then
            yi = sqrt(ei) * cb
        end if 
        if (abs(gi - kci) .gt. ca * gi) then
            kci = 2 * sqrt(ei)
            l = 2 * l
            if (yi .lt. 0) then
                l = l + 1
            end if 
            goto 1
        end if
        if (yi .lt. 0) then 
            l = l + 1
        end if 
        ei = (atan(mi/yi) + pi * l) * (ai / mi)
        if (xi .lt. 0) then
            e(i) = -ei
        else
            e(i) = ei + ci * zi
        end if
2   enddo
    
    return
end

subroutine el3_tmp(x, kc, p, e, j) bind(C, name="el3")

    integer (c_int), bind(C) :: j
    real (c_double), bind(C) :: x(j), kc(j), p(j)
    real (c_double), bind(C), dimension(j), intent(out) :: e
    
    real*8 :: am(j), ap(j), c(j), d(j), de(j), f(j), fa(j), g(j)
    real*8 :: h(j), hh(j), p1(j), pm(j), pz(j), q(j), r(j), s(j) 
    real*8 :: t(j), u(j), v(j), w(j), y(j), ye(j), z(j), zd(j)
    real*8 :: ami, api, ci, di, dei, ei, fi, fai, gi, hi, hhi 
    real*8 :: p1i, pmi, pzi, qi, ri, si, ti, ui, vi, wi, yi, yei, zi, zdi, ppi, xi, kci
    real*8 :: ca, cb
    real*8 :: pi = 3.1415926535897932384626433832795
    real*8 :: ln2 = 0.6931471805599453
    integer :: l, m, n, ND, i, k
    logical :: bo(j)
    logical :: bk, boi
    real*8 :: ra(12), rb(12), rr(12)
    
    ca = 1E-6
    cb = 1E-14
    ND = 10
    
    hh = x * x
    f = p * hh
    
    do i=1,j,1
    
        ami = am(i)
        api = ap(i)
        ci = c(i)
        di = d(i)
        dei = de(i)
        ei = e(i)
        fi = f(i)
        fai = fa(i)
        gi = g(i)
        hi = h(i)
        hhi = hh(i)
        p1i = p1(i)
        pmi = pm(i)
        pzi = pz(i)
        ppi = p(i)
        qi = q(i)
        ri = r(i)
        si = s(i)
        ti = t(i)
        ui = u(i)
        vi = v(i)
        wi = w(i)
        xi = x(i)
        yi = y(i)
        yei = ye(i)
        zi = z(i)
        zdi = zd(i)
        boi = bo(i)
        kci = kc(i)
        
        if (kci == 0.d0) then
            si = ca / (1.d0 + abs(xi))
        else
            si = kci
        end if
        ti = si * si
        pmi = ti * 0.5
        ei = hhi * ti
        zi = abs(fi)
        ri = abs(ppi)
        hi = 1.d0 + hhi
        if ((ei .lt. 1.d0) .AND. (zi .lt. 0.1) .AND. (ti .lt. 1.d0) .AND. (ri .lt. 1.d0)) then
            goto 1
        end if 
        wi = 1.d0 + fi
        if (wi == 0.d0) then
            goto 4
        end if
        if (ppi == 0.d0) then
            p1i = cb / hhi
        else
            p1i = ppi
        end if
        si = abs(si)
        yi = abs(xi)
        gi = p1i - 1.d0
        if (gi == 0.d0) then
            gi = cb
        end if
        fi = p1i - ti
        if (fi == 0.d0) then
            fi = cb * ti
        end if 
        ami = 1.d0 - ti
        api = 1.d0 + ei
        ri = p1i * hi
        fai = gi / (fi * p1i)
        boi = fai .gt. 0.d0
        fai = abs(fai)
        pzi = abs(gi * fi)
        dei = sqrt(pzi)
        qi = sqrt(abs(p1i))
        if (pmi .gt. 0.5) then 
            pmi = 0.5
        end if
        pmi = p1i - pmi
        if (pmi .ge. 0.d0) then
            ui = sqrt(ri * api)
            vi = yi * dei
            if (gi .lt. 0.d0) then
                vi = -vi
            end if
            di = 1.d0 / qi
            ci = 1.d0
        else
            ui = sqrt(hi * api * pzi)
            yei = yi * qi
            vi = ami * yei
            qi = -dei / gi
            di = -ami / dei
            ci = 0.d0
            pzi = api - ri
        end if
        if (boi) then
            ri = vi / ui
            zi = 1.d0
            k = 1
            if (pmi .lt. 0.d0) then
                hi = yi * sqrt(hi / (api * fai))
                hi = 1.d0 / hi - hi
                zi = hi - ri - ri
                ri = 2.d0 + ri * hi
                if (ri == 0.d0) then
                    ri = cb
                end if
                if (zi == 0.d0) then 
                    zi = hi * cb
                end if
                ri = ri / zi
                zi = ri
                wi = pzi
            end if
            ui = ui / wi
            vi = vi / wi
        else
            ti = ui + abs(vi)
            bk = .TRUE.
            if (p1i .lt. 0.d0) then
                dei = vi / pzi
                yei = ui * yei
                yei = yei + yei
                ui = ti / pzi
                vi = (-fi - gi * ei) / ti
                ti = pzi * abs(wi)
                zi = (hhi * ri * fi - gi * api + yei) / ti
                yei = yei / ti
            else
                dei = vi / wi 
                yei = 0.d0
                ui = (ei + p1i) / ti
                vi = ti / wi
                zi = 1.d0
            end if
            if (si .gt. 1.d0) then
                hi = ui
                ui = vi
                vi = hi
            end if
        end if
        yi = 1.d0 / yi
        ei = si
        n = 1
        ti = 1.d0
        m = 0
        l = 0
3       yi = yi - ei / yi
        if (yi == 0.d0) then
            yi = sqrt(ei) * cb
        end if
        fi = ci
        ci = di / qi + ci
        gi = ei / qi
        di = fi * gi + di
        di = di + di
        qi = gi + qi
        gi = ti
        ti = si + ti
        n = n + n
        m = m + m
        if (boi) then 
            if (zi .lt. 0.d0) then 
                m = k + m
            end if 
            k = int(sign(1.d0, ri))
            hi = ei / (ui * ui + vi * vi)
            ui = ui * (1.d0 + hi)
            vi = vi * (1.d0 - hi)
        else
            ri = ui / vi 
            hi = zi * ri
            zi = hi * zi
            hhi = ei / vi
            if (bk) then
                dei = dei / ui
                yei = yei * (hi + 1.d0 / hi) + dei * (1.d0 + ri)
                dei = dei * (ui - hhi)
                bk = abs(yei) .lt. 1.d0
            else
                k = exponent(zi)
                zi = fraction(zi)
                m = m + k
            end if
        end if
        if (abs(gi - si) .gt. ca * gi) then 
            if (boi) then
                gi = (1.d0 / ri - ri) * 0.5
                hhi = ui + vi * gi
                hi = gi * ui - vi
                if (hhi == 0.d0) then
                    hhi = ui * cb
                end if
                if (hi == 0) then
                    hi = vi * cb
                end if
                zi = ri * hi
                ri = hhi / hi
            else
                ui = ui + ei / ui
                vi = vi + hhi
            end if
            si = sqrt(ei)
            si = si + si
            ei = si * ti
            l = l + l
            if (yi .lt. 0.d0) then 
                l = 1 + l
            end if 
            goto 3
        end if
        if (yi .lt. 0.d0) then
            l = 1 + l
        end if
        ei = atan(ti / yi) + pi * l
        ei = ei * (ci * ti + di) / (ti * (ti + qi))
        if (boi) then 
            hi = vi / (ti + ui)
            zi = 1.d0 - ri * hi
            hi = ri + hi
            if (zi == 0.d0) then
                zi = cb
            end if
            if (zi .lt. 0.d0) then
                m = m + int(sign(1.d0, hi))
            end if
            si = atan(hi / zi) + m * pi
        else
            if (bk) then
                si = asin(yei)
            else
                si = log(zi) + m * ln2
            end if
            si = si * 0.5
        end if
        ei = (ei + sqrt(fai) * si) / n
        if (xi .gt. 0.d0) then
            e(i) = ei
        else
            e(i) = -ei
        end if
        goto 4

1       do k=2,ND,1
            rb(k) = 0.5 / k
            ra(k) = 1.d0 - rb(k)
        end do
        zdi = 0.5 / (ND + 1)
        si = ppi + pmi
        do k=2,ND,1
            rr(k) = si
            pmi = pmi * ti * ra(k)
            si = si * ppi + pmi
        end do
        ui = si * zdi
        si = ui
        boi = .FALSE.
        do k=ND,2,-1
            ui = ui + (rr(k) - ui) * rb(k)
            boi = .NOT. boi
            if (boi) then
                vi = -ui
            else
                vi = ui
            end if
            si = si * hhi + vi
        end do
        if (boi) then
            si = -si
        end if
        ui = (ui + 1.d0) * 0.5
        e(i) = (ui - si * hi) * sqrt(hi) * xi + ui * asin(xi)
4   end do
    return
end 

end module bulirsch