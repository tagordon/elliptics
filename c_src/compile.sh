gfortran -c ellip.f90
gfortran -c phot.f90
gfortran -c bulirsch.f90
gcc -O2 -c -Wall -Werror -fpic photlib.c
gfortran ellip.o phot.o bulirsch.o photlib.o -o cwrapper.so