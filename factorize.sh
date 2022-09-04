#!/bin/bash
#shell script for the binary 'lll_factor', which factors a given polynomial input
#passes intput to python script which parses and passes polynomial to c function

argstr="${@:2}"

./bin/lll_factor $(python3 ./utils/poly_parse.py $1) $argstr


