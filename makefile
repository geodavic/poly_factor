all: algebraicmake polymake

algebraicmake:
	gcc -Wall -Wextra -o mpz_algebraic src/mpz_algebraic.c -lgmp -lmpfr -lmpc

polymake:
	gcc -Wall -Wextra -o factor_poly src/factor_poly.c -lgmp -lmpfr -lmpc
