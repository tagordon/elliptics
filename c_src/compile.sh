rm *.o
rm *.mod
gfortran -Ofast -c ellip.f90
gfortran -Ofast -c phot.f90
gfortran -Ofast -c kepler.f90
gcc -O3 -c -Wall -Werror -fpic photlib.c
gfortran -Ofast ellip.o phot.o kepler.o photlib.o -o cwrapper.so
