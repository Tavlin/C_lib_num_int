#include "initial_data_and_functions.h"

//left Riemann sum
double left_riemann_sum(InitialData A, Gaussian params,
double(*func)(double, Gaussian));


//right Rieman sum
double right_riemann_sum(InitialData A, Gaussian params,
double(*func)(double, Gaussian));

//Trapezodial Rule
double trapezodial_integral(InitialData A, Gaussian params,
double(*func)(double, Gaussian));
