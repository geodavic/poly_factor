# Polynomial factorization using LLL

This is a repository implementing the LLL lattice reduction algorithm and using it for monic, integer polynomial factorization. Given such a polynomial p, the algorithm works roughly like this: use a quickly converging root finding method (e.g. cubic Newton's method) to find a root p(a) = 0, where a is possibly complex. Then try to find the minimal polynomial p_a of a, which must be an irreducible factor of p. This is done by finding an integer relation on the powers of a using the lattice basis reduction algorithm due to Lenstra, Lenstra, and Lovaz (LLL). Divide p by p_a and repeat this process on the quotient q = p/p_a. Eventually the process will stop when q = 1. 

The higher the degree of the polynomial, the higher the precision needed; hence the GMP dependency. 

## Use

Use this [page](https://web.ma.utexas.edu/users/gdavtor/poly_factor.html) to try it out. 

If you prefer the command line, clone this repo and use the included makefile to get two binaries: ```factor_poly``` and ```mpz_algebraic```. In other words:
```
make all //compile both files and run unit tests
```

The first binary is polynomial factorization; for example:
```
bin/lll_factor -1,0,0,0,1
```
will factor the polynomial x^4-1. 

The second is a algebraic number checker. Execute the binary and follow the inputs. It will take as input a base 10 decimal number and output a feasible minimal polynomial for that number. 

There is also a shell script version of ```lll_factor``` located in the root directory that can be used in the terminal for easier input of the polynomial, e.g.:
```
./factorize.sh x^4-1
```
This will have the same output as ```bin/lll_factor -1,0,0,0,1```. 

## Dependencies

Compiling requires the libraries ```gmp 5.+```, ```mpfr 1.1+```, and ```mpc 3.+```. 

## Unit tests

The following is a good unit test for development. It consists of factoring a large degree polynomial that has repeated factors.
```
./factorize.sh x^18+3x^17+6x^16+2x^15+5x^14+9x^13-5x^12-3x^10-17x^9-3x^7-16x^6+5x^5+9x^4-5x^3+3x^2+6x -newline

Polynomial input: 6x + 3x^2 - 5x^3 + 9x^4 + 5x^5 - 16x^6 - 3x^7 - 17x^9 - 3x^10 - 5x^12 + 9x^13 + 5x^14 + 2x^15 + 6x^16 + 3x^17 + x^18
Factorization:

(-1 - x + x^3 + x^4 + x^5)
(1 - x + x^2)^3
(1 + x)^3
(-1 + x)
(6 + 3x + x^2)
(x)
```

## API

There is an API hosted on Google Cloud which you can use to factor polynomials using your favorite requests engine. For example, with cURL:
```
curl -d '{"poly":"x^4+x^2+1"}' -H "Content-Type: application/json" -X POST "https://poly-factor-s4ph7avbaq-uc.a.run.app/factor"
```
The full documentation for this API can be found [here](https://poly-factor-s4ph7avbaq-uc.a.run.app/docs)

## Dockerfile

There is a Dockerfile included that can be used to compile the code should you not have the libraries. Simply build the image and exec into it, and the binaries will be built and ready to use.

## TODO

1) Optimize using ```pthreads``` and generalize to non-monic factorization.

