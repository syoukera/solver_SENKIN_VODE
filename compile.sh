#!/bin/bash

gfortran -o src/dmath.o   -c src/dmath.f -fPIC -fno-range-check
gfortran -o src/twopnt.o  -c src/twopnt.f
gfortran -o src/tranlib.o -c src/tranlib.f
gfortran -o src/cklib.o   -c src/cklib.f
gfortran -c main.f90

gfortran -o solve_senkin main.o src/cklib.o src/tranlib.o src/dmath.o