#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "integrals.h"

//planed:
//changing all given parameters to pointers. Gonna happen... someday... I guess
// printf function to give it a better look and fix some double posting

// if lower boundary > upper boundary
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

// function to print solutions of Integrals
void printer(const char * func_name, const char * integral_name,
double integral_val)
{
	printf("\n*************************************************************\n");
	printf("Integraltype:\t%s\nFunction:\t%s\nIntegralvalue = %+6.10lf\n",
	integral_name, func_name, integral_val);
	printf("*************************************************************\n");
}


// function trafo for open boundary
/*
double * (double,  FunctionParams)func_trafo(double(*func)(double, FunctionParams))
{	
	double x;
	FunctionParams params;
	return 1/pow(x,2) * (*func)(1/x, params);
}*/

//left Riemann sum
double left_riemann_sum(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps)
//initial instead of init, cuz __init__ is SPECIAL, at least somewhere else :)

{
	int swaped;
	// checks that initial value is allways smaller then final value 
	if(A.final_val < A.initial)
	{	
		swaped = 1;
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
	
	if(swaped ==1)
	{
		return (-left_sum);
	}
	
	else
	{
		return left_sum;
	}
}


//right Rieman sum
double right_riemann_sum(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps)
{

	int swaped;
	// checks that initial value is allways smaller then final value 
	if(A.final_val < A.initial)
	{	
		swaped = 1;
		swap(A.final_val, A.initial);
	}
	
	double right_sum = 0;
	
	double stepsize = (A.final_val-A.initial)/A.N;
	
	for(double i = 1; i <= A.N; i++)
	{
		right_sum += stepsize *
		(*func)(((A.initial+(i*(A.final_val-A.initial))/A.N)), params);
	}
	
	if(swaped ==1)
	{
		return (-right_sum);
	}
	
	else
	{
		return right_sum;
	}
}

//Trapezodial Rule
double trapezodial_integral(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps)
{

	int swaped;
	// checks that initial value is allways smaller then final value 
	if(A.final_val < A.initial)
	{
		swaped = 1;
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
	
	if(swaped ==1)
	{
		return (-integral_val);
	}
	
	else
	{
		return integral_val;
	}
	
}

//Simpson's Rule
double simpson_integral(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps)
{

	int swaped;
	// checks that initial value is allways smaller then final value 
	if(A.final_val < A.initial)
	{
		swaped = 1;
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
	
	if(swaped ==1)
	{
		return (-integral_val);
	}
	
	else
	{
		return integral_val;
	}

}


//Trapezodial Rule semi adaptive stepsize
double trapezodial_integral_sas(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps)
{
	int swaped;
	// checks that initial value is allways smaller then final value 
	if(A.final_val < A.initial)
	{
		swaped = 1;
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
	
	printf("Number of steps used in last calculation = %.0lf\n\n", A.N);
	
	if(swaped ==1)
	{
		return (-integral_val);
	}
	
	else
	{
		return integral_val;
	}
	
	
}


// semi adaptive midpoint rule
double midpoint_int(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps)
{
	int swaped;
	// checks that initial value is allways smaller then final value 
	if(A.final_val < A.initial)
	{
		swaped = 1;
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
	
	printf("Number of steps used in last calculation = %.0lf\n\n", A.N);
	
	if(swaped ==1)
	{
		return (-integral_val);
	}
	
	else
	{
		return integral_val;
	}


	
	
}

// semi adaptive midpoint rule with [a,inf)
double midpoint_int_to_inf(InitialData  A, FunctionParams params,
double(*func)(double, FunctionParams), double eps)
{
	
	// check if lower bound is < 0, so split is needed in [a,0.25] nad [0.25,inf)
	// 0.25 so the calculation for both integrals does not take too long/too many
	// steps
	
	double f_int;
	
	if(A.initial <= 0)
	{
		A.final_val = 0.25;
		f_int = midpoint_int(A, params,(*func), eps);
		A.initial = 0;
		A.final_val = 4;
	}
	
	if(A.initial > 0)
	{
		A.final_val = 1/A.initial;
		A.initial = 0;
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
	(*p_previous_int_val) = h/pow(midpoint,2) * (*func)(1.0/midpoint, params);
	
	// so that N is up to date with N = 3 BEFORE entering the while loop!
	A.N *= 3.0;
	h /= 3.0;
	
	// initial int value N = 3
	(*p_integral_val) = (*p_previous_int_val)/3.0 +  h/pow(midpoint-h,2)*
	((*func)(midpoint - h, params)) + h/pow(midpoint+h,2) *
	((*func)(midpoint + h, params));
	
	while(fabs((*p_previous_int_val)-(*p_integral_val)) > eps)
	{
		
		A.N *= 3.0;
		h /= 3.0;
	
		(*p_previous_int_val) = (*p_integral_val);
	
		*p_new_int_val = 0;
	
		for(double j = -((A.N/2.0) - 0.5); j <= ((A.N/2.0) - 0.5); j += 3)
		{
			(*p_new_int_val) += h/pow((midpoint + (j*h)),2) * 
			(*func)(1/(midpoint + (j*h)), params);
			
		}
	
		for(double i = (-((A.N/2.0) - 0.5)) + 2; i <= ((A.N/2.0) - 0.5); i += 3)
		{
			(*p_new_int_val) += h/pow((midpoint + (i*h)),2) * 
			(*func)(1/(midpoint + (i*h)), params);
			
		}
	
		(*p_integral_val) = (*p_new_int_val) + (*p_integral_val)/3.0;

	}
	integral_val += f_int;
	printf("Number of steps used in last calculation = %.0lf\n\n", A.N);
	
	return integral_val;
	
}


