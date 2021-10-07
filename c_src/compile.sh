rm *.o
rm *.mod
gfortran -Ofast -unroll=12 -c ellip.f90
gfortran -Ofast -c -unroll=5 phot.f90
gfortran -Ofast -g -shared ellip.o phot.o -o fwrapper.so
gfortran -Ofast -c -unroll=5 kepler.f90
gcc -O3 -c -Wall -Werror -fpic photlib.c
gfortran -Ofast ellip.o phot.o kepler.o photlib.o -o cwrapper.so
