all: algebraic poly unit

algebraic:
	gcc -Wall -Wextra -o mpz_algebraic src/mpz_algebraic.c -lgmp -lmpfr -lmpc

poly:
	gcc -Wall -Wextra -o factor_poly src/factor_poly.c -lgmp -lmpfr -lmpc

unit:
	src/factorize.sh `cat test_poly.txt`
