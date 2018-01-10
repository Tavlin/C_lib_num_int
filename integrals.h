#include "initial_data_and_functions.h"

// Swap algorithm
void xorSwap (double x, double y);

// function to print solutions of Integrals
void printer(const char * func_name, const char * integral_name,
double integral_val);

// function trafo for open boundary
//double func_trafo(double(*func)(double x, FunctionParams));


// left Riemann sum
double left_riemann_sum(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps);

// right Riemann sum
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

// semi adaptive midpoint rule with [a,inf)
double midpoint_int_to_inf(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps);
