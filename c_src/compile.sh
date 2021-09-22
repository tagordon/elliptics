rm *.o
rm *.mod
gfortran -O3 -c ellip.f90
gfortran -O3 -c phot.f90
gfortran -O3 -c kepler.f90
gcc -O3 -c -Wall -Werror -fpic photlib.c
gfortran -O3 ellip.o phot.o kepler.o photlib.o -o cwrapper.so
