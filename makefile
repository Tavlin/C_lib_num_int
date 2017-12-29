prog: main.c initial_data.c integrals.c
	gcc -o prog main.c initial_data.c integrals.c -I. -lm
