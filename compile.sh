#!/bin/bash

FORTRAN=gfortran

FLAGS="-I/usr/include  -Wall -fopenmp  -O3 -march=native -ffpe-trap=invalid,overflow,zero"
DFLAGS="-I/usr/include  -Wall -fopenmp -ffpe-trap=invalid,overflow,zero -fcheck=bounds -fbacktrace -g"

rm -f srelax.x 

#                                         FFTW3     LAPACK
 $FORTRAN  $FLAGS srelax.f90 -o srelax.x -lfftw3 # -llapack -lblas
#$FORTRAN $DFLAGS srelax.f90 -o srelax.x -lfftw3 # -llapack -lblas
