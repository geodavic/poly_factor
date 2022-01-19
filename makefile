build: algebraic poly
units: unit1 unit2
all: algebraic poly unit1 unit2

algebraic:
	gcc -Wall -Wextra -o bin/mpz_algebraic src/mpz_algebraic.c -lgmp -lmpfr -lmpc

poly:
	gcc -Wall -Wextra -o bin/factor_poly src/factor_poly.c -lgmp -lmpfr -lmpc

unit1:
	./factorize.sh `cat test/test_poly.txt` > test/make_output.txt
	diff test/make_output.txt test/correct_output.txt

unit2:
	./factorize.sh `cat test/test_poly2.txt` > test/make_output2.txt
	diff test/make_output2.txt test/correct_output2.txt
