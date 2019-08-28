#!/bin/bash
#shell script for the binary 'factorize', which factors a given polynomial input
#passes intput to python script which parses and passes polynomial to c function
#five possible arguments for now

argstr="${@:2}"

rootstr=/Users/David/Documents/Maths/c/algebraic_numbers/factorization
$rootstr/commandline_factor/factor_poly $(python $rootstr/msc_code/poly_parse.py $1) $argstr


