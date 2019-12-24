all: algebraic poly unit1 unit2
units: unit1 unit2

algebraic:
	gcc -Wall -Wextra -o bin/mpz_algebraic src/mpz_algebraic.c -lgmp -lmpfr -lmpc

poly:
	gcc -Wall -Wextra -o bin/factor_poly src/factor_poly.c -lgmp -lmpfr -lmpc

unit1:
	./factorize.sh `cat unit_tests/test_poly.txt` > unit_tests/make_output.txt
	diff unit_tests/make_output.txt unit_tests/correct_output.txt

unit2:
	./factorize.sh `cat unit_tests/test_poly2.txt` > unit_tests/make_output2.txt
	diff unit_tests/make_output2.txt unit_tests/correct_output2.txt
