prog: main.c initial_data_and_functions.c integrals.c 
	gcc -o prog main.c initial_data_and_functions.c integrals.c -I. -lm
