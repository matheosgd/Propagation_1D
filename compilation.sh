#!/bin/sh

FLAGS="-Og -g -fbacktrace -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -finit-real=nan -finit-integer=nan -fopenmp"

MODEXT="-IQDUtilLib-main/OBJ/obj_gfortran_opt0_omp1_lapack1"
LIB="QDUtilLib-main/libQD_gfortran_opt0_omp1_lapack1.a   -llapack -lblas"

rm Main.exe
gfortran -c Basis.f90 $FLAGS $MODEXT
gfortran -c Operateurs.f90 $FLAGS $MODEXT
gfortran -c Psi.f90 $FLAGS $MODEXT
gfortran -c Lanczos.f90 $FLAGS $MODEXT
gfortran -c Propagation.f90 $FLAGS $MODEXT
gfortran -c Main.f90 $FLAGS $MODEXT
gfortran -o Main.exe $FLAGS Basis.o Operateurs.o Psi.o Lanczos.o Propagation.o TD_Matheo.o  $LIB

#gfortran -Og -g -fbacktrace -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -finit-real=nan -finit-integer=nan TD_Matheo.f90 -o TD_Matheo.exe

./Main.exe < data_morse > results