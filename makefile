build-c: algebraic-c poly-c
build: build-c

test-factor: unit1 unit2
test: test-factor test-poly-cc

all: build test

algebraic-c:
	gcc -Wall -Wextra -o bin/mpz_algebraic src/c/mpz_algebraic.c -lgmp -lmpfr -lmpc

poly-c:
	gcc -Wall -Wextra -o bin/lll_factor src/c/lll_factor.c -lgmp -lmpfr -lmpc

unit1:
	@echo "----------------------- Factor test 1 ------------------------"
	@./factorize.sh `cat test/factor/test_poly.txt` > test/factor/make_output.txt
	diff test/factor/make_output.txt test/factor/correct_output.txt

unit2:
	@echo "----------------------- Factor test 2 ------------------------"
	@./factorize.sh `cat test/factor/test_poly2.txt` > test/factor/make_output2.txt
	diff test/factor/make_output2.txt test/factor/correct_output2.txt

test-poly-cc:
	@g++ -std=c++17 -I ./src/cpp test/cpp/test_poly.cc -lgmpxx -lgmp -o tester.o
	@./tester.o
	@rm tester.o
