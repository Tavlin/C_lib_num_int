prog: const_stepsize_alg.c initial_data.c integrals.c
	gcc -o prog const_stepsize_alg.c initial_data.c integrals.c -I. -lm
