gcc -Wall -c heat_eqn.c
gcc heat_eqn.o -lgsl -lgslcblas -lm
./a.out
