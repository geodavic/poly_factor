# Polynomial factorization using LLL

This is a repository implementing the LLL lattice reduction algorithm and using it for monic, integer polynomial factorization. 

## Use

Use the included makefile to make two binaries: ```factor_poly``` and ```mpz_algebraic```. The first is the polynomial factorization; for example:
```
./factor_poly -1,0,0,0,1
```
will factor the polynomial x^4-1. The second is a algebraic number checker; to run it, just execute the binary and follow the inputs. It will take as input a base 10 decimal number and output a feasible minimal polynomial for that number. 
