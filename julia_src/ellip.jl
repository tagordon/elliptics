 #module elli2
#implicit none

#!real*8 :: el1, el2, el3

#contains

function el1(x::T, kc::T) where {T <: Real}

#    integer :: j
#    real*8 :: x, kc, e
#    !real, intent(out) :: e
    
#    real*8 :: g, m, y
#    real*8 :: D = 12
    D = 12
#    real*8 :: ca, cb
#    real*8 :: pi = 3.1415926535897932384626433832795
#    integer :: l, i
    
    ca = exp10(-D/2.0)
    cb = exp10(-D+2.0)
    
    y = abs(1.0/x)
    kc = abs(kc)
    m = 1.0
    
    if (kc == 0.0)
        ee = log(x + 1.0 / cos(atan(x)))  # 1/cos(atan(x)) = sqrt(1+x^2)
        @goto two
    end
    
    l = 0
        
@label one
    ee = m * kc
    g = m
    m = kc + m
    y = - (ee/y) + y
    if (y == 0)
        y = sqrt(ee) * cb
    end
    if (abs(g - kc) > (ca * g))
        kc = 2 * sqrt(ee)
        l = 2 * l
        if (y < 0)
            l += 1
        end
        @goto one
    end
    if (y < 0)
        l += 1
    end
    ee = (atan(m,y) + pi * l) / m
    if (x < 0)
        ee = -ee
    end
    
@label two
   el1 = ee
return el1
end

function el2(x::T, kc::T, a::T, b::T) where {T <: Real}

#    real*8 :: x, kc, a, b, e
#    !real, intent(out) :: e
    
#    real*8 :: c, dd, f, g, ik, m, p, y, z
#    real*8 :: D = 12
    D = 12
#    real*8 :: ca, cb
#    real*8 :: pi = 3.1415926535897932384626433832795
#    integer :: l, i
    
    ca = exp10(-D/2.0)
    cb = exp10(-D+2.0)
    
    c = x^2
    dd = c + 1.0
    p = sqrt((1.0 + kc^2 * c)/dd)
    dd = x / dd
    c = dd / (2.0 * p)
    z = a - b
    ik = a
    a = (b + a) * 0.5
    y = abs(1.0/x)
    f = 0.0
    kc = abs(kc)
    m = 1.0
        
    if (kc == 0.0)
        ee = sin(atan(x))  # Isn't this x/sqrt(1+x^2)?
        @goto two 
    end

    l = 0
@label one
    b = ik * kc + b
    ee = m * kc
    g = ee / p
    dd = f * g + dd
    f = c
    ik = a
    p = g + p
    c = (dd / p + c) / 2.0
    g = m
    m = kc + m
    a = (b / m + a) * 0.5
    y = - (ee/y) + y
    if (y == 0)
        y = sqrt(ee) * cb
    end
    if (abs(g - kc) > ca * g)
        kc = 2 * sqrt(ee)
        l = 2 * l
        if (y < 0)
            l = l + 1
        end
        @goto one
    end
    if (y < 0)
        l = l + 1
    end
    ee = (atan(m/y) + pi * l) * (a / m)
    if (x < 0)
        ee = -ee
    else
        ee = ee + c * z
    end
    
@label two
    el2 = ee
return el2
end

function el3(x::T, kc::T, p::T) where {T <: Real}

#    real*8 :: x, kc, p, e
#    !real, intent(out) :: e
    
#    real*8 :: am, ap, c, d, de, f, fa, g
#    real*8 :: h, hh, p1, pm, pz, q, r, s 
#    real*8 :: t, u, v, w, y, ye, z, zd
#    real*8 :: ca, cb
#    real*8 :: pi = 3.1415926535897932384626433832795
#    real*8 :: ln2 = 0.6931471805599453
    ln2 = log(convert(T,2))
#    integer :: l, m, n, ND, i, k
#    logical :: bo, bk
#    real*8 :: ra(10), rb(10), rr(10)
    ra=zeros(10); rb=zeros(10); rr=zeros(10)
    
    ca = 1e-6
    cb = 1e-14
    ND = 10      # Should these be higher for higher precision?
    
    hh = x * x
    f = p * hh
        
    if (kc == 0.0)
        s = ca / (1.0 + abs(x))
    else
        s = kc
    end
    t = s * s
    pm = t * 0.5
    ee = hh * t
    z = abs(f)
    r = abs(p)
    h = 1.0 + hh
    if ((ee < 1.0) && (z < 0.1) && (t < 1.0) && (r < 1.0))
        @goto one
    end
    w = 1.0 + f
    if (w == 0.0)
        @goto four
    end
    if (p == 0.0)
        p1 = cb / hh
    else
        p1 = p
    end
    s = abs(s)
    y = abs(x)
    g = p1 - 1.0
    if (g == 0.0)
        g = cb
    end
    f = p1 - t
    if (f == 0.0)
        f = cb * t
    end
    am = 1.0 - t
    ap = 1.0 + ee
    r = p1 * h
    fa = g / (f * p1)
    bo = fa > 0.0
    fa = abs(fa)
    pz = abs(g * f)
    de = sqrt(pz)
    q = sqrt(abs(p1))
    if (pm > 0.5)
        pm = 0.5
    end
    pm = p1 - pm
    if (pm >= 0.0)
        u = sqrt(r * ap)
        v = y * de
        if (g < 0.0)
            v = -v
        end
        d = 1.0 / q
        c = 1.0
    else
        u = sqrt(h * ap * pz)
        ye = y * q
        v = am * ye
        q = -de / g
        d = -am / de
        c = 0.0
        pz = ap - r
    end
    if (bo)
        r = v / u
        z = 1.0
        k = 1
        if (pm < 0.0)
            h = y * sqrt(h / (ap * fa))
            h = 1.0 / h - h
            z = h - r - r
            r = 2.0 + r * h
            if (r == 0.0)
                r = cb
            end
            if (z == 0.0)
                z = h * cb
            end
            r = r / z
            z = r
            w = pz
        end
        u = u / w
        v = v / w
    else
        t = u + abs(v)
        bk = true
        if (p1 < 0.0)
            de = v / pz
            ye = u * ye
            ye = ye + ye
            u = t / pz
            v = (-f - g * ee) / t
            t = pz * abs(w)
            z = (hh * r * f - g * ap + ye) / t
            ye = ye / t
        else
            de = v / w 
            ye = 0.0
            u = (ee + p1) / t
            v = t / w
            z = 1.0
        end
        if (s > 1.0)
            h = u
            u = v
            v = h
        end
    end
    y = 1.0 / y
    ee = s
    n = 1
    t = 1.0
    m = 0
    l = 0
@label three
    y = y - ee / y
    if (y == 0.0)
        y = sqrt(ee) * cb
    end
    f = c
    c = d / q + c
    g = ee / q
    d = f * g + d
    d = d + d
    q = g + q
    g = t
    t = s + t
    n = n + n
    m = m + m
    if (bo)
        if (z < 0.0)
            m += k
        end
        k = convert(Int64,sign(r))  # Check syntax
        h = ee / (u * u + v * v)
        u = u * (1.0 + h)
        v = v * (1.0 - h)
    else
        r = u / v 
        h = z * r
        z = h * z
        hh = ee / v
        if (bk)
            de = de / u
            ye = ye * (h + 1.0 / h) + de * (1.0 + r)
            de = de * (u - hh)
            bk = abs(ye) < 1.0
        else
            k = exponent(z)+1  # Check syntax
            z = fraction(z)/2  # Check syntax
            m = m + k
        end
    end
    if (abs(g - s) > ca * g)
        if (bo)
            g = (1.0 / r - r) * 0.5
            hh = u + v * g
            h = g * u - v
            if (hh == 0.0)
                hh = u * cb
            end
            if (h == 0.0)
                h = v * cb
            end
            z = r * h
            r = hh / h
        else
            u = u + ee / u
            v = v + hh
        end
        s = sqrt(ee)
        s = s + s
        ee = s * t
        l = l + l
        if (y < 0.0)
            l += 1
        end
        @goto three 
    end
    if (y < 0.0)
        l += 1 
    end
    ee = atan(t,y) + pi * l   # Maybe atan(t,y)?
    ee = ee * (c * t + d) / (t * (t + q))
    if (bo)
        h = v / (t + u)
        z = 1.0 - r * h
        h = r + h
        if (z == 0.0)
            z = cb
        end
        if (z < 0.0)
            m = m + convert(Int64,sign(h))  # Check syntax
        end
        s = atan(h,z) + m * pi
    else
        if (bk)
            s = asinh(ye)
        else
            s = log(z) + m * ln2
        end
        s = s * 0.5
    end
    ee = (ee + sqrt(fa) * s) / n
    if (x <= 0.0)
        ee = -ee
    end
    @goto four 

@label one
    for k=2:ND 
        rb[k] = 0.5 / k
        ra[k] = 1.0 - rb[k]
    end
    zd = 0.5 / (ND + 1)
    s = p + pm
    for k=2:ND
        rr[k] = s
        pm = pm * t * ra[k]
        s = s * p + pm
    end
    u = s * zd
    s = u
    bo = false
    for k=ND:-1:2
        u = u + (rr[k] - u) * rb[k]
        bo = !bo
        if (bo)
            v = -u
        else
            v = u
        end
        s = s * hh + v
    end
    if (bo)
        s = -s
    end
    u = (u + 1.0) * 0.5
    ee = (u - s * h) * sqrt(h) * x + u * asinh(x)
    
@label four
    el3 = ee
return el3
end

# This routine can be found in the Limbdark.jl repository
# as cel_bulirsch
function cel(kc::T, p::T, a::T, b::T) where {T <: Real}
#
#    real*8 :: kc, p, a, b
#    real*8 :: CA, e, f, g, m, q
#    
  CA = 1e-6
    
  kc = abs(kc)
  ee = kc
  m = 1.0
    
  if kc == 0.0
    cel = 0.0 
    return cel
  end
    
  if p > 0.0
    p = sqrt(p)
    b = b / p
  else
    f = kc * kc
    q = 1.0 - f
    g = 1.0 - p
    f = f - p
    q = (b - a * p) * q
    p = sqrt(f / g)
    a = (a - b) / g
    b = -q / (g * g * p) + a * p
  end

@label one
  f = a
  a = b / p + a
  g = ee / p
  b = f * g + b
  b = b + b
  p = g + p
  g = m
  m = kc + m
    
  if abs(g - kc) > (g * CA)
    kc = sqrt(ee)
    kc = kc + kc
    ee = kc * m
    @goto one
  end
  cel = 0.5*pi* ((a * m + b) / (m * (m + p)))

return cel
end

#end module ellip
