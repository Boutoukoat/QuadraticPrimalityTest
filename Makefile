

GGG = g++ -O3 -Wall -march=native -fomit-frame-pointer -fexpensive-optimizations

OBJ = quadratic_primality_main.o \
      quadratic_primality.o \
      quadratic_primality_alloc.o \
      expression_parser.a

quadratic: $(OBJ)
	$(GGG) -static -o quadratic $(OBJ) -lgmp -lpthread

quadratic_primality_main.o: quadratic_primality_main.cpp quadratic_primality.h quadratic_primality_alloc.h expression_parser.h
	$(GGG) -c -o quadratic_primality_main.o quadratic_primality_main.cpp

quadratic_primality_alloc.o: quadratic_primality_alloc.cpp quadratic_primality_alloc.h
	$(GGG) -c -o quadratic_primality_alloc.o quadratic_primality_alloc.cpp

quadratic_primality.o: quadratic_primality.cpp quadratic_primality.h
	$(GGG) -c -o quadratic_primality.o quadratic_primality.cpp

expression_parser.a : bison.gmp_expr.o lex.gmp_expr.o expression_parser.h
	ar vr expression_parser.a bison.gmp_expr.o lex.gmp_expr.o

bison.gmp_expr.o : bison.gmp_expr.tab.c bison.gmp_expr.h
	$(GGG) -c -o bison.gmp_expr.o bison.gmp_expr.tab.c

bison.gmp_expr.tab.c bison.gmp_expr.tab.h : parser.y
	bison -d parser.y

lex.gmp_expr.o : lex.gmp_expr.c
	$(GGG) -Wno-unused-function -c -o lex.gmp_expr.o lex.gmp_expr.c

lex.gmp_expr.c : parser.l bison.gmp_expr.tab.h
	flex parser.l

check: quadratic
	./quadratic -st

clean:
	rm -f ./quadratic $(OBJ)


