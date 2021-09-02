rm *.o
rm *.mod
gfortran -c ellip.f90
gfortran -c phot.f90
gcc -O2 -c -Wall -Werror -fpic photlib.c
gfortran ellip.o phot.o photlib.o -o cwrapper.so
