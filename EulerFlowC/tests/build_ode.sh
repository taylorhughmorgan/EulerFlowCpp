gcc -Wall -c ode_example.c
gcc ode_example.o -lgsl -lgslcblas -lm
./a.out
