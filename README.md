# Polynomial factorization using LLL

This is a repository implementing the LLL lattice reduction algorithm and using it for monic, integer polynomial factorization. 

## Use

Use the included makefile to compile two binaries: ```factor_poly``` and ```mpz_algebraic```. The first is polynomial factorization; for example:
```
./factor_poly -1,0,0,0,1
```
will factor the polynomial x^4-1. The second is a algebraic number checker. Execute the binary and follow the inputs. It will take as input a base 10 decimal number and output a feasible minimal polynomial for that number. 

There is also a shell script version of ```factor_poly``` located in ```src``` that can be used in the terminal for easier input of the polynomial, e.g.:
```
factorize.sh x^4-1
```
This will have the same output as ```./factor_poly -1,0,0,0,1```. 

## Dependencies

Compiling requires the libraries ```gmp 5.+```, ```mpfr 1.1+```, and ```mpc 3.+```. 
