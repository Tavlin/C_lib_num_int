#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "integrals.h"

//planed:
//changing all given parameters to pointers. Gonna happen... someday... I guess

// https://en.wikipedia.org/wiki/XOR_swap_algorithm
void swap (double x, double y)
{
	if (x != y)
	{
		double b;
		x = b;
		x = y;
		y = b;
	}
}


//left Riemann sum
double left_riemann_sum(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps)
//initial instead of init, cuz __init__ is SPECIAL, at least somewhere else :)

{
	// checks that initial value is allways smaller then final value 
	if(A.final_val < A.initial)
	{
		swap(A.final_val, A.initial);
	}
	
	// start value for the sum
	double left_sum = 0;
	
	double stepsize = (A.final_val-A.initial)/A.N; 
	
	for(double i = 1; i <= A.N; i++)
	{
		left_sum += stepsize *
		(*func)((A.initial + (((i-1)*(A.final_val-A.initial))/A.N)), params);
	}
	printf("left riemann sum = %+6.10lf\n", left_sum);
}


//right Rieman sum
double right_riemann_sum(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps)
{

	// checks that initial value is allways smaller then final value 
	if(A.final_val < A.initial)
	{
		swap(A.final_val, A.initial);
	}
	
	double right_sum = 0;
	
	double stepsize = (A.final_val-A.initial)/A.N;
	
	for(double i = 1; i <= A.N; i++)
	{
		right_sum += stepsize *
		(*func)(((A.initial+(i*(A.final_val-A.initial))/A.N)), params);
	}
	printf("right riemann sum = %+6.10lf\n", right_sum);
}

//Trapezodial Rule
double trapezodial_integral(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps)
{

	// checks that initial value is allways smaller then final value 
	if(A.final_val < A.initial)
	{
		swap(A.final_val, A.initial);
	}
	
	// upper function value with pointy pointer
	double right_step = 0;
	double * pright_step;
	pright_step = &right_step;
	
	// lower function value with pointy pointer
	double left_step = 0;
	double * pleft_step;
	pleft_step = &left_step;
	
	double integral_val = 0;
	
	double stepsize = (A.final_val-A.initial)/A.N;
	
	*pleft_step = (*func)((A.initial), params);
	
	for(double i = 1; i <= A.N; i++)
	{
		*pright_step = (*func)((A.initial+(i*(A.final_val-A.initial))/A.N), params);
		
		integral_val += stepsize/2.0 * 
		((*pleft_step) + (*pright_step));
		
		(*pleft_step) = (*pright_step);
	}
	printf("trapezodial integral = %+6.10lf\n", integral_val);
}

//Simpson's Rule
double simpson_integral(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps)
{

	// checks that initial value is allways smaller then final value 
	if(A.final_val < A.initial)
	{
		swap(A.final_val, A.initial);
	}
	
	// upper function value with pointy pointer
	double right_step = 0;
	double * pright_step;
	pright_step = &right_step;
	
	// lower function value with pointy pointer
	double left_step = 0;
	double * pleft_step;
	pleft_step = &left_step;
	
	double integral_val = 0;
	
	double stepsize = (A.final_val-A.initial)/A.N;
	
	// initial lower function value
	*pleft_step = (*func)((A.initial), params);
	
	for(double i = 1; i <= A.N; i++)
	{
		*pright_step = (*func)((A.initial+(i*(A.final_val-A.initial))/A.N), params);
		
		integral_val += (stepsize/6.0) *
		((*pleft_step) + (*pright_step) + 
		(4.0*(*func)(((A.initial+(i*(A.final_val-A.initial))/A.N) +
		(A.initial + (((i-1)*(A.final_val-A.initial))/A.N)))/2.0,params)));
		
		(*pleft_step) = (*pright_step);
		
	}
	printf("simpson integral = %+6.10lf\n\n", integral_val);
}


//Trapezodial Rule semi adaptive stepsize
double trapezodial_integral_sas(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps)
{

	// checks that initial value is allways smaller then final value 
	if(A.final_val < A.initial)
	{
		swap(A.final_val, A.initial);
	}
	
	// reset # steps to 1
	A.N = 1.0;
	
	double integral_val = 0;
	double* p_integral_val = &integral_val;
	
	double previous_int_val = 0;
	double* p_previous_int_val = &previous_int_val;
	
	double stepsize = 1.0/A.N;
	double h = (A.final_val - A.initial);
	
	double new_int_val = 0;
	double* p_new_int_val = &new_int_val;
	
	// initial int value N = 1
	(*p_previous_int_val) = (h/2.0) *
	((*func)(A.initial,params) + (*func)(A.final_val,params));

	// so that N is up to date with N = 2 BEFORE entering the while loop!
	A.N *= 2.0;
	h /= 2.0;

	// initial int value N = 2
	(*p_integral_val) = (*p_previous_int_val)/2.0 + 
	h * ((*func)((A.final_val + A.initial)/2.0, params));
	

	while(fabs((*p_previous_int_val)-(*p_integral_val)) > eps)
	{
		
		A.N *= 2.0;
		h /= 2.0;
		stepsize = 1.0/A.N;
	
		// 2*stepsize so the calculation is only done once per stepsize
		// otherwise it would need to calc for every for iteration
		double tt_stepsize = 2.0*stepsize;
	
		// new maximum value for the new midpoints:
		double max_val_midpoints = (A.N-1.0)/A.N;
	
		(*p_previous_int_val) = (*p_integral_val);
	
		*p_new_int_val = 0;
	
		for(double j = 1; j <= A.N; j += 2)
		{
			(*p_new_int_val) += h * 
			(*func)((A.initial + j*h), params);
			
		}
	
		(*p_integral_val) = (*p_new_int_val) + (*p_integral_val)/2.0;

	}
	printf("Trapezodial Rule with semi adaptive stepsize = %+6.10lf\n",
	integral_val);
	printf("Number of steps used in last calculation = %.0lf\n\n", A.N);
	
}


// semi adaptive midpoint rule
double midpoint_int(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps)
{

	// checks that initial value is allways smaller then final value 
	if(A.final_val < A.initial)
	{
		swap(A.final_val, A.initial);
	}
	
	// reset # steps to 1
	A.N = 1.0;
	
	double integral_val = 0;
	double* p_integral_val = &integral_val;
	
	double previous_int_val = 0;
	double* p_previous_int_val = &previous_int_val;
	
	double stepsize = 1.0/A.N;
	double h = (A.final_val - A.initial);
	
	double midpoint = (A.final_val + A.initial)/2.0;
	
	double new_int_val = 0;
	double* p_new_int_val = &new_int_val;
	
	// initial int value N = 1
	(*p_previous_int_val) = h * (*func)(midpoint, params);
	
	// so that N is up to date with N = 3 BEFORE entering the while loop!
	A.N *= 3.0;
	h /= 3.0;
	
	// initial int value N = 3
	(*p_integral_val) = (*p_previous_int_val)/3.0 + h * 
	(((*func)(midpoint - h, params)) + ((*func)(midpoint + h, params)));
	
	//printf("f(-h) = %+6.10lf\n", (*func)(midpoint - h, params));
	//printf("f(+h) = %+6.10lf\n", (*func)(midpoint + h, params));
	//printf("Integral for N = 1 is %+6.10lf\n", (*p_previous_int_val));
	//printf("Integral for N = 3 is %+6.10lf\n", (*p_integral_val));
	
	while(fabs((*p_previous_int_val)-(*p_integral_val)) > eps)
	{
		
		A.N *= 3.0;
		h /= 3.0;
	
		(*p_previous_int_val) = (*p_integral_val);
	
		*p_new_int_val = 0;
	
		for(double j = -((A.N/2.0) - 0.5); j <= ((A.N/2.0) - 0.5); j += 3)
		{
			(*p_new_int_val) += h * (*func)(midpoint + (j*h), params);
			
		}
	
		for(double i = (-((A.N/2.0) - 0.5)) + 2; i <= ((A.N/2.0) - 0.5); i += 3)
		{
			(*p_new_int_val) += h * (*func)(midpoint + (i*h), params);
			
		}
	
		(*p_integral_val) = (*p_new_int_val) + (*p_integral_val)/3.0;

	}
	printf("Midpoint Rule with semi adaptive stepsize = %+6.10lf\n",
	integral_val);
	printf("Number of steps used in last calculation = %.0lf\n\n", A.N);
	
	
}






















