#!/bin/bash
#shell script for the binary 'factor_poly', which factors a given polynomial input
#passes intput to python script which parses and passes polynomial to c function

argstr="${@:2}"
rootstr="/Users/David/Documents/Maths/c/algebraic_numbers/git"

$rootstr/factor_poly $(python $rootstr/msc_code/poly_parse.py $1) $argstr


