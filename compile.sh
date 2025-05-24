#!/bin/bash

FORTRAN=gfortran

FLAGS="-Wall -fopenmp"

rm -f qrelax.x 

$FORTRAN $FLAGS -o qrelax.x  qrelax.f90
