module ellip
use iso_c_binding
implicit none

contains

subroutine FKP_C(phi, n, m, f, e, p, k) bind(C, name="fkp")

    real*8 :: pi = 3.14159265358979323856264338327950288419716939937510
    integer :: i
    integer (C_INT), BIND(C) :: k
    real (C_DOUBLE), BIND(C) :: phi(k), n(k), m(k)
    real*8 :: mc(k), s(k), ncirc(k)
    real (C_DOUBLE), BIND(C), dimension(k), intent(out) :: f, e, p
    real*8 :: btmp, dtmp, jtmp, mc_inv, n_inv, s_inv
    
    mc = 1.d0 - m
    s = sin(phi)

    do i=1,k,1
        
        if (m(i) .lt. 1.d0) then
            call elsbdj(s(i), n(i), mc(i), btmp, dtmp, jtmp)    
            f(i) = (btmp + dtmp)
            e(i) = (btmp + mc(i)*dtmp)
            p(i) = (btmp + dtmp + n(i)*jtmp)
        else
            mc_inv = 1.d0 - 1.d0/m(i)
            n_inv = n(i) / m(i)
            s_inv = s(i) / sqrt(m(i))
            call elsbdj(s_inv, n_inv, mc_inv, btmp, dtmp, jtmp)
            f(i) = (btmp + dtmp)
            e(i) = (btmp + mc(i)*dtmp)
            p(i) = (btmp + dtmp + n(i)*jtmp)
        end if
    enddo

    return 
end

subroutine elsbdj_arrays(s0, n, m, b, d, j, l)

    integer :: i, flag
    integer :: l
    real*8 :: s0(l), n(l), m(l), mc(l)
    real*8, dimension(l), intent(out) :: b, d, j
    real*8 btmp, dtmp, jtmp, sqm
    
    mc = 1.d0 - m

    do i=1,l,1
        
        call elsbdj(s0(i), n(i), mc(i), btmp, dtmp, jtmp)
        
        b(i) = btmp
        d(i) = dtmp
        j(i) = jtmp
    enddo

    return 
end

subroutine elsbdj_arrays_c(s0, n, m, b, d, j, k) bind(C, name="bdj_arrays")

integer :: i, flag
integer (C_INT), bind(C) :: k
real (C_DOUBLE), bind(C) :: s0(k), n(k), m(k)
real*8 ::  mc(k)
real (C_DOUBLE), dimension(k), intent(out) :: b, d, j
real*8 btmp, dtmp, jtmp, sqm, nm, nn, ns0
mc = 1.d0 - m

do i=1,k,1
    call elsbdj(s0(i), n(i), mc(i), btmp, dtmp, jtmp)
    
    b(i) = btmp
    d(i) = dtmp
    j(i) = jtmp
enddo

return 
end


subroutine elsbdj(s0,n,mc,b,d,j)
!
! Simultaneous computation of associate elliptic integrals
! of third kind, B(phi|m), D(phi|m), and J(phi,n|m)
! by using the half/double argument transformation of sn functions
!
! Reference: Fukushima, T, (2012) J. Comp. Appl. Math., 236, 1961-1975
!   Precise and fast computation of a general incomplete elliptic integral
!   of third kind by half and double argument transformations
!
real*8 s0, n, mc
real*8 b, d, j
real*8 m, h, del, s, y, c, sy, t
real*8 yy(11),ss(11),cd(11)

!for calling from c, change this to: real (C_DOUBLE), bind(C) :: serj, uatan
!real*8 serj, uatan
integer i,k

m=1.d0-mc
h=n*(1.d0-n)*(n-m)
del=0.01622d0 

s=s0
y=s*s
if(y.lt.del) then
    call serbd(y,m,b,d)
    b=s*b
    d=s*y*d
    j=s*serj(y,n,m)
    return
endif
yy(1)=y
ss(1)=s
do i=1,10
    c=sqrt(1.d0-y)
    d=sqrt(1.d0-m*y)
    y=y/((1.d0+c)*(1.d0+d))
    yy(i+1)=y
    ss(i+1)=sqrt(y)
    cd(i)=c*d
    if(y.lt.del) then
        goto 1
    endif
enddo

1 continue
call serbd(y,m,b,d)
b=ss(i+1)*b
d=ss(i+1)*y*d
j=ss(i+1)*serj(y,n,m)
do k=i,1,-1
    sy=ss(k)*yy(k+1)
    t=sy/(1.d0-n*(yy(k)-yy(k+1)*cd(k)))
    b=2.d0*b-sy
    d=d+(d+sy)
    j=j+(j+uatan(t,h))
enddo
return
end

!---------------------------------------------------------------------------
subroutine serbd(y,m,b,d)
!
! Simultaneous computation of associate elliptic integrals,
! B(phi|m) and D(phi|m), for small arguments by the series expansion
!
! Reference: Fukushima, T (2012) J. Comp. Appl. Math., 235, 4140-4148
!   Precise and fast computation of general incomplete elliptic integral
!   of second kind by half and double argument transformations
!
real*8 y,m,b,d
real*8 F1,F2,F3,F4
real*8 F10,F20,F21,F30,F31,F40,F41,F42
real*8 F5,F50,F51,F52,F6,F60,F61,F62,F63
real*8 F7,F70,F71,F72,F73,F8,F80,F81,F82,F83,F84
real*8 F9,F90,F91,F92,F93,F94
real*8 FA,FA0,FA1,FA2,FA3,FA4,FA5
real*8 FB,FB0,FB1,FB2,FB3,FB4,FB5
parameter (F10=1.d0/6.d0)
parameter (F20=3.d0/40.d0)
parameter (F21=2.d0/40.d0)
parameter (F30=5.d0/112.d0)
parameter (F31=3.d0/112.d0)
parameter (F40=35.d0/1152.d0)
parameter (F41=20.d0/1152.d0)
parameter (F42=18.d0/1152.d0)
parameter (F50=63.d0/2816.d0)
parameter (F51=35.d0/2816.d0)
parameter (F52=30.d0/2816.d0)
parameter (F60=231.d0/13312.d0)
parameter (F61=126.d0/13312.d0)
parameter (F62=105.d0/13312.d0)
parameter (F63=100.d0/13312.d0)
parameter (F70=429.d0/30720.d0)
parameter (F71=231.d0/30720.d0)
parameter (F72=189.d0/30720.d0)
parameter (F73=175.d0/30720.d0)
parameter (F80=6435.d0/557056.d0)
parameter (F81=3432.d0/557056.d0)
parameter (F82=2722.d0/557056.d0)
parameter (F83=2520.d0/557056.d0)
parameter (F84=2450.d0/557056.d0)
parameter (F90=12155.d0/1245184.d0)
parameter (F91=6435.d0/1245184.d0)
parameter (F92=5148.d0/1245184.d0)
parameter (F93=4620.d0/1245184.d0)
parameter (F94=4410.d0/1245184.d0)
parameter (FA0=46189.d0/5505024.d0)
parameter (FA1=24310.d0/5505024.d0)
parameter (FA2=19305.d0/5505024.d0)
parameter (FA3=17160.d0/5505024.d0)
parameter (FA4=16170.d0/5505024.d0)
parameter (FA5=15876.d0/5505024.d0)
parameter (FB0=88179.d0/12058624.d0)
parameter (FB1=46189.d0/12058624.d0)
parameter (FB2=36465.d0/12058624.d0)
parameter (FB3=32175.d0/12058624.d0)
parameter (FB4=30030.d0/12058624.d0)
parameter (FB5=29106.d0/12058624.d0)
real*8 A1,A2,A3,A4,A5,A6,A7,A8,A9,AA,AB
parameter (A1=3.d0/5.d0)
parameter (A2=5.d0/7.d0)
parameter (A3=7.d0/9.d0)
parameter (A4=9.d0/11.d0)
parameter (A5=11.d0/13.d0)
parameter (A6=13.d0/15.d0)
parameter (A7=15.d0/17.d0)
parameter (A8=17.d0/19.d0)
parameter (A9=19.d0/21.d0)
parameter (AA=21.d0/23.d0)
parameter (AB=23.d0/25.d0)
real*8 B1,B2,B3,B4,B5,B6,B7,B8,B9,BA,BB
real*8 D0,D1,D2,D3,D4,D5,D6,D7,D8,D9,DA,DB
parameter (D0=1.d0/3.d0)
!
!	write(*,*) "(serbd) y,m=",y,m
F1=F10+m*F10
F2=F20+m*(F21+m*F20)
F3=F30+m*(F31+m*(F31+m*F30))
F4=F40+m*(F41+m*(F42+m*(F41+m*F40)))
F5=F50+m*(F51+m*(F52+m*(F52+m*(F51+m*F50))))
F6=F60+m*(F61+m*(F62+m*(F63+m*(F62+m*(F61+m*F60)))))
F7=F70+m*(F71+m*(F72+m*(F73+m*(F73+m*(F72+m*(F71+m*F70))))))
F8=F80+m*(F81+m*(F82+m*(F83+m*(F84+m*(F83+m*(F82+m*(F81+m*F80)))))))
F9=F90+m*(F91+m*(F92+m*(F93+m*(F94+m*(F94+m*(F93 &
    +m*(F92+m*(F91+m*F90))))))))
FA=FA0+m*(FA1+m*(FA2+m*(FA3+m*(FA4+m*(FA5+m*(FA4 &
    +m*(FA3+m*(FA2+m*(FA1+m*FA0)))))))))
FB=FB0+m*(FB1+m*(FB2+m*(FB3+m*(FB4+m*(FB5+m*(FB5+ &
    m*(FB4+m*(FB3+m*(FB2+m*(FB1+m*FB0))))))))))
!
D1=F1*A1
D2=F2*A2
D3=F3*A3
D4=F4*A4
D5=F5*A5
D6=F6*A6
D7=F7*A7
D8=F8*A8
D9=F9*A9
DA=FA*AA
DB=FB*AB
d=D0+y*(D1+y*(D2+y*(D3+y*(D4+y*(D5+y*(D6+y*(D7+y*(D8 &
    +y*(D9+y*(DA+y*DB))))))))))
B1=F1-D0
B2=F2-D1
B3=F3-D2
B4=F4-D3
B5=F5-D4
B6=F6-D5
B7=F7-D6
B8=F8-D7
B9=F9-D8
BA=FA-D9
BB=FB-DA
b=1.d0+y*(B1+y*(B2+y*(B3+y*(B4+y*(B5+y*(B6+y*(B7+y*(B8 &
    +y*(B9+y*(BA+y*BB))))))))))
return
end
!---------------------------------------------------------------------------
! add bind(C) statement to call from c
real*8 function serj(y,n,m)
!
! Computation of associate elliptic integral J(phi,n|m)
! for small arguments by the series expansion
!
! Reference: Fukushima, T, (2012) J. Comp. Appl. Math., 236, 1961-1975
!   Precise and fast computation of a general incomplete elliptic integral
!   of third kind by half and double argument transformations
!
real*8 y,n,m

real*8 J1,J2,J3,J4,J5,J6,J7,J8,J9,JA

real*8 J100,J200,J201,J210,J300,J301,J302,J310,J311,J320
real*8 J400,J401,J402,J403,J410,J411,J412,J420,J421,J430
real*8 J500,J501,J502,J503,J504,J510,J511,J512,J513,J520
real*8 J521,J522,J530,J531,J540
real*8 J600,J601,J602,J603,J604,J605,J610,J611,J612,J613,J614
real*8 J620,J621,J622,J623,J630,J631,J632,J640,J641,J650
real*8 J700,J701,J702,J703,J704,J705,J706
real*8 J710,J711,J712,J713,J714,J715,J720,J721,J722,J723,J724
real*8 J730,J731,J732,J733,J740,J741,J742,J750,J751,J760
real*8 J800,J801,J802,J803,J804,J805,J806,J807
real*8 J810,J811,J812,J813,J814,J815,J816
real*8 J820,J821,J822,J823,J824,J825,J830,J831,J832,J833,J834
real*8 J840,J841,J842,J843,J850,J851,J852,J860,J861,J870
real*8 J900,J901,J902,J903,J904,J905,J906,J907,J908
real*8 J910,J911,J912,J913,J914,J915,J916,J917
real*8 J920,J921,J922,J923,J924,J925,J926
real*8 J930,J931,J932,J933,J934,J935,J940,J941,J942,J943,J944
real*8 J950,J951,J952,J953,J960,J961,J962,J970,J971,J980
real*8 JA00,JA01,JA02,JA03,JA04,JA05,JA06,JA07,JA08,JA09
real*8 JA10,JA11,JA12,JA13,JA14,JA15,JA16,JA17,JA18
real*8 JA20,JA21,JA22,JA23,JA24,JA25,JA26,JA27
real*8 JA30,JA31,JA32,JA33,JA34,JA35,JA36
real*8 JA40,JA41,JA42,JA43,JA44,JA45,JA50,JA51,JA52,JA53,JA54
real*8 JA60,JA61,JA62,JA63,JA70,JA71,JA72,JA80,JA81,JA90

parameter (J100=1.d0/3.d0)

parameter (J200=1.d0/10.d0)
parameter (J201=2.d0/10.d0)
parameter (J210=1.d0/10.d0)

parameter (J300=3.d0/56.d0)
parameter (J301=4.d0/56.d0)
parameter (J302=8.d0/56.d0)
parameter (J310=2.d0/56.d0)
parameter (J311=4.d0/56.d0)
parameter (J320=3.d0/56.d0)

parameter (J400=5.d0/144.d0)
parameter (J401=6.d0/144.d0)
parameter (J402=8.d0/144.d0)
parameter (J403=16.d0/144.d0)
parameter (J410=3.d0/144.d0)
parameter (J411=4.d0/144.d0)
parameter (J412=8.d0/144.d0)
parameter (J420=3.d0/144.d0)
parameter (J421=6.d0/144.d0)
parameter (J430=5.d0/144.d0)

parameter (J500=35.d0/1408.d0)
parameter (J501=40.d0/1408.d0)
parameter (J502=48.d0/1408.d0)
parameter (J503=64.d0/1408.d0)
parameter (J504=128.d0/1408.d0)
parameter (J510=20.d0/1408.d0)
parameter (J511=24.d0/1408.d0)
parameter (J512=32.d0/1408.d0)
parameter (J513=64.d0/1408.d0)
parameter (J520=18.d0/1408.d0)
parameter (J521=24.d0/1408.d0)
parameter (J522=48.d0/1408.d0)
parameter (J530=20.d0/1408.d0)
parameter (J531=40.d0/1408.d0)
parameter (J540=35.d0/1408.d0)

parameter (J600=63.d0/3328.d0)
parameter (J601=70.d0/3328.d0)
parameter (J602=80.d0/3328.d0)
parameter (J603=96.d0/3328.d0)
parameter (J604=128.d0/3328.d0)
parameter (J605=256.d0/3328.d0)
parameter (J610=35.d0/3328.d0)
parameter (J611=40.d0/3328.d0)
parameter (J612=48.d0/3328.d0)
parameter (J613=64.d0/3328.d0)
parameter (J614=128.d0/3328.d0)
parameter (J620=30.d0/3328.d0)
parameter (J621=36.d0/3328.d0)
parameter (J622=48.d0/3328.d0)
parameter (J623=96.d0/3328.d0)
parameter (J630=30.d0/3328.d0)
parameter (J631=40.d0/3328.d0)
parameter (J632=80.d0/3328.d0)
parameter (J640=35.d0/3328.d0)
parameter (J641=70.d0/3328.d0)
parameter (J650=63.d0/3328.d0)

parameter (J700=231.d0/15360.d0)
parameter (J701=252.d0/15360.d0)
parameter (J702=280.d0/15360.d0)
parameter (J703=320.d0/15360.d0)
parameter (J704=384.d0/15360.d0)
parameter (J705=512.d0/15360.d0)
parameter (J706=1024.d0/15360.d0)
parameter (J710=126.d0/15360.d0)
parameter (J711=140.d0/15360.d0)
parameter (J712=160.d0/15360.d0)
parameter (J713=192.d0/15360.d0)
parameter (J714=256.d0/15360.d0)
parameter (J715=512.d0/15360.d0)
parameter (J720=105.d0/15360.d0)
parameter (J721=120.d0/15360.d0)
parameter (J722=144.d0/15360.d0)
parameter (J723=192.d0/15360.d0)
parameter (J724=384.d0/15360.d0)
parameter (J730=100.d0/15360.d0)
parameter (J731=120.d0/15360.d0)
parameter (J732=160.d0/15360.d0)
parameter (J733=320.d0/15360.d0)
parameter (J740=105.d0/15360.d0)
parameter (J741=140.d0/15360.d0)
parameter (J742=280.d0/15360.d0)
parameter (J750=126.d0/15360.d0)
parameter (J751=252.d0/15360.d0)
parameter (J760=231.d0/15360.d0)

parameter (J800=429.d0/34816.d0)
parameter (J801=462.d0/34816.d0)
parameter (J802=504.d0/34816.d0)
parameter (J803=560.d0/34816.d0)
parameter (J804=640.d0/34816.d0)
parameter (J805=768.d0/34816.d0)
parameter (J806=1024.d0/34816.d0)
parameter (J807=2048.d0/34816.d0)
parameter (J810=231.d0/34816.d0)
parameter (J811=252.d0/34816.d0)
parameter (J812=280.d0/34816.d0)
parameter (J813=320.d0/34816.d0)
parameter (J814=384.d0/34816.d0)
parameter (J815=512.d0/34816.d0)
parameter (J816=1024.d0/34816.d0)
parameter (J820=189.d0/34816.d0)
parameter (J821=210.d0/34816.d0)
parameter (J822=240.d0/34816.d0)
parameter (J823=288.d0/34816.d0)
parameter (J824=284.d0/34816.d0)
parameter (J825=768.d0/34816.d0)
parameter (J830=175.d0/34816.d0)
parameter (J831=200.d0/34816.d0)
parameter (J832=240.d0/34816.d0)
parameter (J833=320.d0/34816.d0)
parameter (J834=640.d0/34816.d0)
parameter (J840=175.d0/34816.d0)
parameter (J841=210.d0/34816.d0)
parameter (J842=280.d0/34816.d0)
parameter (J843=560.d0/34816.d0)
parameter (J850=189.d0/34816.d0)
parameter (J851=252.d0/34816.d0)
parameter (J852=504.d0/34816.d0)
parameter (J860=231.d0/34816.d0)
parameter (J861=462.d0/34816.d0)
parameter (J870=429.d0/34816.d0)

parameter (J900=6435.d0/622592.d0)
parameter (J901=6864.d0/622592.d0)
parameter (J902=7392.d0/622592.d0)
parameter (J903=8064.d0/622592.d0)
parameter (J904=8960.d0/622592.d0)
parameter (J905=10240.d0/622592.d0)
parameter (J906=12288.d0/622592.d0)
parameter (J907=16384.d0/622592.d0)
parameter (J908=32768.d0/622592.d0)
parameter (J910=3432.d0/622592.d0)
parameter (J911=3696.d0/622592.d0)
parameter (J912=4032.d0/622592.d0)
parameter (J913=4480.d0/622592.d0)
parameter (J914=5120.d0/622592.d0)
parameter (J915=6144.d0/622592.d0)
parameter (J916=8192.d0/622592.d0)
parameter (J917=16384.d0/622592.d0)
parameter (J920=2772.d0/622592.d0)
parameter (J921=3024.d0/622592.d0)
parameter (J922=3360.d0/622592.d0)
parameter (J923=3840.d0/622592.d0)
parameter (J924=4608.d0/622592.d0)
parameter (J925=6144.d0/622592.d0)
parameter (J926=12288.d0/622592.d0)
parameter (J930=2520.d0/622592.d0)
parameter (J931=2800.d0/622592.d0)
parameter (J932=3200.d0/622592.d0)
parameter (J933=3840.d0/622592.d0)
parameter (J934=5120.d0/622592.d0)
parameter (J935=10240.d0/622592.d0)
parameter (J940=2450.d0/622592.d0)
parameter (J941=2800.d0/622592.d0)
parameter (J942=3360.d0/622592.d0)
parameter (J943=4480.d0/622592.d0)
parameter (J944=8960.d0/622592.d0)
parameter (J950=2520.d0/622592.d0)
parameter (J951=3024.d0/622592.d0)
parameter (J952=4032.d0/622592.d0)
parameter (J953=8064.d0/622592.d0)
parameter (J960=2772.d0/622592.d0)
parameter (J961=3696.d0/622592.d0)
parameter (J962=7392.d0/622592.d0)
parameter (J970=3432.d0/622592.d0)
parameter (J971=6864.d0/622592.d0)
parameter (J980=6435.d0/622592.d0)

parameter (JA00=12155.d0/1376256.d0)
parameter (JA01=12870.d0/1376256.d0)
parameter (JA02=13728.d0/1376256.d0)
parameter (JA03=14784.d0/1376256.d0)
parameter (JA04=16128.d0/1376256.d0)
parameter (JA05=17920.d0/1376256.d0)
parameter (JA06=20480.d0/1376256.d0)
parameter (JA07=24576.d0/1376256.d0)
parameter (JA08=32768.d0/1376256.d0)
parameter (JA09=65536.d0/1376256.d0)
parameter (JA10=6435.d0/1376256.d0)
parameter (JA11=6864.d0/1376256.d0)
parameter (JA12=7392.d0/1376256.d0)
parameter (JA13=8064.d0/1376256.d0)
parameter (JA14=8960.d0/1376256.d0)
parameter (JA15=10240.d0/1376256.d0)
parameter (JA16=12288.d0/1376256.d0)
parameter (JA17=16384.d0/1376256.d0)
parameter (JA18=32768.d0/1376256.d0)
parameter (JA20=5148.d0/1376256.d0)
parameter (JA21=5544.d0/1376256.d0)
parameter (JA22=6048.d0/1376256.d0)
parameter (JA23=6720.d0/1376256.d0)
parameter (JA24=7680.d0/1376256.d0)
parameter (JA25=9216.d0/1376256.d0)
parameter (JA26=12288.d0/1376256.d0)
parameter (JA27=24576.d0/1376256.d0)
parameter (JA30=4620.d0/1376256.d0)
parameter (JA31=5040.d0/1376256.d0)
parameter (JA32=5600.d0/1376256.d0)
parameter (JA33=6400.d0/1376256.d0)
parameter (JA34=7680.d0/1376256.d0)
parameter (JA35=10240.d0/1376256.d0)
parameter (JA36=20480.d0/1376256.d0)
parameter (JA40=4410.d0/1376256.d0)
parameter (JA41=4900.d0/1376256.d0)
parameter (JA42=5600.d0/1376256.d0)
parameter (JA43=6720.d0/1376256.d0)
parameter (JA44=8960.d0/1376256.d0)
parameter (JA45=17920.d0/1376256.d0)
parameter (JA50=4410.d0/1376256.d0)
parameter (JA51=5040.d0/1376256.d0)
parameter (JA52=6048.d0/1376256.d0)
parameter (JA53=8064.d0/1376256.d0)
parameter (JA54=16128.d0/1376256.d0)
parameter (JA60=4620.d0/1376256.d0)
parameter (JA61=5544.d0/1376256.d0)
parameter (JA62=7392.d0/1376256.d0)
parameter (JA63=14784.d0/1376256.d0)
parameter (JA70=5148.d0/1376256.d0)
parameter (JA71=6864.d0/1376256.d0)
parameter (JA72=13728.d0/1376256.d0)
parameter (JA80=6435.d0/1376256.d0)
parameter (JA81=12870.d0/1376256.d0)
parameter (JA90=12155.d0/1376256.d0)

! write(*,"(a20,1p3e10.2)") "(serj) y,n,m=",y,n,m

J1=J100
J2=J200+n*J201+m*J210

J3=J300+n*(J301+n*J302)+m*(J310+n*J311+m*J320)

J4=J400+n*(J401+n*(J402+n*J403)) &
    +m*(J410+n*(J411+n*J412)+m*(J420+n*J421+m*J430))

J5=J500+n*(J501+n*(J502+n*(J503+n*J504))) &
    +m*(J510+n*(J511+n*(J512+n*J513)) &
    +m*(J520+n*(J521+n*J522)+m*(J530+n*J531+m*J540)))
if(y.le.6.0369310d-04) then
    serj=y*(J1+y*(J2+y*(J3+y*(J4+y*J5))))
    return
endif

J6=J600+n*(J601+n*(J602+n*(J603+n*(J604+n*J605)))) &
    +m*(J610+n*(J611+n*(J612+n*(J613+n*J614))) &
    +m*(J620+n*(J621+n*(J622+n*J623)) &
    +m*(J630+n*(J631+n*J632)+m*(J640+n*J641+m*J650))))
if(y.le.2.0727505d-03) then
    serj=y*(J1+y*(J2+y*(J3+y*(J4+y*(J5+y*J6)))))
! write(*,"(a20,1pe10.2)") "(serj) J6",serj
    return
endif

J7=J700+n*(J701+n*(J702+n*(J703+n*(J704+n*(J705+n*J706))))) &
    +m*(J710+n*(J711+n*(J712+n*(J713+n*(J714+n*J715)))) &
    +m*(J720+n*(J721+n*(J722+n*(J723+n*J724))) &
    +m*(J730+n*(J731+n*(J732+n*J733)) &
    +m*(J740+n*(J741+n*J742)+m*(J750+n*J751+m*J760)))))
if(y.le.5.0047026d-03) then
    serj=y*(J1+y*(J2+y*(J3+y*(J4+y*(J5+y*(J6+y*J7))))))
   return
endif

J8=J800+n*(J801+n*(J802+n*(J803+n*(J804+n*(J805+n*(J806+n*J807)))))) &
    +m*(J810+n*(J811+n*(J812+n*(J813+n*(J814+n*(J815+n*J816))))) &
    +m*(J820+n*(J821+n*(J822+n*(J823+n*(J824+n*J825)))) &
    +m*(J830+n*(J831+n*(J832+n*(J833+n*J834))) &
    +m*(J840+n*(J841+n*(J842+n*J843)) &
    +m*(J850+n*(J851+n*J852)+m*(J860+n*J861+m*J870))))))
if(y.le.9.6961652d-03) then
    serj=y*(J1+y*(J2+y*(J3+y*(J4+y*(J5+y*(J6+y*(J7+y*J8)))))))
    return
endif

J9=J900+n*(J901+n*(J902+n*(J903+n*(J904+n*(J905+n*(J906+n*(J907+n*J908))))))) &
    +m*(J910+n*(J911+n*(J912+n*(J913+n*(J914+n*(J915+n*(J916+n*J917)))))) &
    +m*(J920+n*(J921+n*(J922+n*(J923+n*(J924+n*(J925+n*J926))))) &
    +m*(J930+n*(J931+n*(J932+n*(J933+n*(J934+n*J935)))) &
    +m*(J940+n*(J941+n*(J942+n*(J943+n*J944))) &
    +m*(J950+n*(J951+n*(J952+n*J953)) &
    +m*(J960+n*(J961+n*J962)+m*(J970+n*J971+m*J980)))))))
if(y.le.1.6220210d-02) then
    serj=y*(J1+y*(J2+y*(J3+y*(J4+y*(J5+y*(J6+y*(J7+y*(J8+y*J9))))))))
    return
endif

JA=JA00+n*(JA01+n*(JA02+n*(JA03+n*(JA04+n*(JA05+n*(JA06+n*(JA07+n*(JA08+n*JA09)))))))) &
    +m*(JA10+n*(JA11+n*(JA12+n*(JA13+n*(JA14+n*(JA15+n*(JA16+n*(JA17+n*JA18))))))) &
    +m*(JA20+n*(JA21+n*(JA22+n*(JA23+n*(JA24+n*(JA25+n*(JA26+n*JA27)))))) &
    +m*(JA30+n*(JA31+n*(JA32+n*(JA33+n*(JA34+n*(JA35+n*JA36))))) &
    +m*(JA40+n*(JA41+n*(JA42+n*(JA43+n*(JA44+n*JA45)))) &
    +m*(JA50+n*(JA51+n*(JA52+n*(JA53+n*JA54))) &
    +m*(JA60+n*(JA61+n*(JA62+n*JA63)) &
    +m*(JA70+n*(JA71+n*JA72)+m*(JA80+n*JA81+m*JA90))))))))
    serj=y*(J1+y*(J2+y*(J3+y*(J4+y*(J5+y*(J6+y*(J7+y*(J8+y*(J9+y*JA)))))))))
    return

end
!---------------------------------------------------------------------------
! add bind(C) statement to call from c
real*8 function uatan(t,h)
!
! Universal arctangent function
!
! Reference: Fukushima, T, (2012) J. Comp. Appl. Math., 236, 1961-1975
!   Precise and fast computation of a general incomplete elliptic integral
!   of third kind by half and double argument transformations
!
real*8 t,h,z,y,x
real*8 a,r,ri,hold,rold,riold
real*8 A3,A5,A7,A9,A11,A13,A15,A17,A19,A21,A23,A25
data hold/1.d0/, rold/1.d0/,riold/1.d0/
save hold,rold,riold
parameter (A3=1.d0/3.d0)
parameter (A5=1.d0/5.d0)
parameter (A7=1.d0/7.d0)
parameter (A9=1.d0/9.d0)
parameter (A11=1.d0/11.d0)
parameter (A13=1.d0/13.d0)
parameter (A15=1.d0/15.d0)
parameter (A17=1.d0/17.d0)
parameter (A19=1.d0/19.d0)
parameter (A21=1.d0/21.d0)
parameter (A23=1.d0/23.d0)
parameter (A25=1.d0/25.d0)

z=-h*t*t
a=abs(z)

if(a.lt.3.3306691d-16) then
    uatan=t
elseif(a.lt.2.3560805d-08) then
    uatan=t*(1.d0+z*A3)
elseif(a.lt.9.1939631d-06) then
    uatan=t*(1.d0+z*(A3+z*A5))
elseif(a.lt.1.7779240d-04) then
    uatan=t*(1.d0+z*(A3+z*(A5+z*A7)))
elseif(a.lt.1.0407839d-03) then
    uatan=t*(1.d0+z*(A3+z*(A5+z*(A7+z*A9))))
elseif(a.lt.3.3616998d-03) then
    uatan=t*(1.d0+z*(A3+z*(A5+z*(A7+z*(A9+z*A11)))))
elseif(a.lt.7.7408014d-03) then
    uatan=t*(1.d0+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*A13))))))
elseif(a.lt.1.4437181d-02) then
    uatan=t*(1.d0+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*A15)))))))
elseif(a.lt.2.3407312d-02) then
    uatan=t*(1.d0+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*(A15+z*A17))))))))
elseif(a.lt.3.4416203d-02) then
    uatan=t*(1.d0+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*(A15+z*(A17+z*A19)))))))))
elseif(z.lt.0.d0) then
    if(abs(h-hold).lt.1.d-16) then
        r=rold; ri=riold
    else
        r=sqrt(h); ri=1.d0/r; hold=h; rold=r; riold=ri
    endif
    uatan=atan(r*t)*ri
elseif(a.lt.4.7138547d-02) then
    uatan=t*(1.d0+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*(A15+z*(A17+z*(A19+z*A21))))))))))
elseif(a.lt.6.1227405d-02) then
    uatan=t*(1.d0+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*(A15+z*(A17+z*(A19+z*(A21+z*A23)))))))))))
elseif(a.lt.7.6353468d-02) then
    uatan=t*(1.d0+z*(A3+z*(A5+z*(A7+z*(A9+z*(A11+z*(A13+z*(A15+z*(A17+z*(A19+z*(A21+z*(A23+z*A25))))))))))))
else
    if(abs(h-hold).lt.1.d-16) then
        r=rold; ri=riold
    else
        r=sqrt(-h); ri=1.d0/r; hold=h; rold=r; riold=ri
    endif
    y=r*t
    x=log((1.d0+y)/(1.d0-y))*0.5d0
    uatan=x*ri
endif

return
end

end module ellip