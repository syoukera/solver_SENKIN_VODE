#!/bin/bash

gfortran -o src/dmath.o   -c src/dmath.f -fPIC -fno-range-check
gfortran -o src/twopnt.o  -c src/twopnt.f
gfortran -o src/tranlib.o -c src/tranlib.f
gfortran -o src/cklib.o   -c src/cklib.f
gfortran -o src/chemkin_m.o -c src/chemkin_m.f90 -J src
gfortran -o src/main.o -c main.f90 -I src

gfortran -o solve_senkin src/main.o src/chemkin_m.o src/cklib.o src/tranlib.o src/dmath.o