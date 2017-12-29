#include "initial_data_and_functions.h"

//left Riemann sum
double left_riemann_sum(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams));


//right Rieman sum
double right_riemann_sum(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams));

//Trapezodial Rule
double trapezodial_integral(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams));

//Simpson's Rule
double simpson_integral(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams));
