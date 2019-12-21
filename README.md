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

## Unit test

The following is a good unit test for development. It consists of factoring a large degree polynomial that has repeated factors.
```
./factorize.sh x^18+3x^17+6x^16+2x^15+5x^14+9x^13-5x^12-3x^10-17x^9-3x^7-16x^6+5x^5+9x^4-5x^3+3x^2+6x

Polynomial input:
6x^1 + 3x^2 - 5x^3 + 9x^4 + 5x^5 - 16x^6 - 3x^7 - 17x^9 - 3x^10 - 5x^12 + 9x^13 + 5x^14 + 2x^15 + 6x^16 + 3x^17 + x^18
Factorization:
x
(-1 - x^1 + x^3 + x^4 + x^5)
(1 - x^1 + x^2)
(1 + x^1)
(-1 + x^1)
(6 + 3x^1 + x^2)
(1 - x^1 + x^2)
(1 + x^1)
(1 - x^1 + x^2)
(1 + x^1)
```
