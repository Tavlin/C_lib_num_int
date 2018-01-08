#include "initial_data_and_functions.h"

// xorSwap algorithm
void xorSwap (double x, double y);

// left Riemann sum
double left_riemann_sum(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps);


// right Rieman sum
double right_riemann_sum(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps);

// Trapezodial Rule
double trapezodial_integral(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps);

// Simpson's Rule
double simpson_integral(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps);

// Trapezodial Rule semi adaptive stepsize
double trapezodial_integral_sas(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps);

// semi adaptive midpoint rule
double midpoint_int(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps);
