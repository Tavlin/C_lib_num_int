prog: const_stepsize_alg.c initial_data.c
	gcc -o prog const_stepsize_alg.c initial_data.c -I. -lm
