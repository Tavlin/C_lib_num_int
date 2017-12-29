#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "integrals.h"

//left Riemann sum
double left_riemann_sum(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams))
//initial instead of init, cuz __init__ is SPECIAL, at least somewhere else :)

{
	double left_sum = 0;
	
	double stepsize;
	stepsize = (A.final_val-A.initial)/A.N;
	
	for(double i = 1; i <= A.N; i++)
	{
		left_sum += stepsize *
		(*func)((A.initial + (((i-1)*(A.final_val-A.initial))/A.N)), params);
	}
	printf("left riemann sum = %+6.10lf\n", left_sum);
}


//right Rieman sum
double right_riemann_sum(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams))
{
	double right_sum = 0;
	
	double stepsize;
	stepsize = (A.final_val-A.initial)/A.N;
	
	for(double i = 1; i <= A.N; i++)
	{
		right_sum += stepsize *
		(*func)(((A.initial+(i*(A.final_val-A.initial))/A.N)), params);
	}
	printf("right riemann sum = %+6.10lf\n", right_sum);
}

//Trapezodial Rule
double trapezodial_integral(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams))
{
	double right_step = 0;
	double * pright_step;
	pright_step = &right_step;
	
	double left_step = 0;
	double * pleft_step;
	pleft_step = &left_step;
	
	double integral_val = 0;
	double stepsize;
	stepsize = (A.final_val-A.initial)/A.N;
	*pleft_step = (*func)((A.initial), params);
	
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
double(*func)(double, FunctionParams))
{
	double right_step = 0;
	double * pright_step;
	pright_step = &right_step;
	
	double left_step = 0;
	double * pleft_step;
	pleft_step = &left_step;
	
	double integral_val = 0;
	double stepsize;
	stepsize = (A.final_val-A.initial)/A.N;
	*pleft_step = (*func)((A.initial), params);
	
	for(double i = 1; i <= A.N; i++)
	{
		*pright_step = (*func)((A.initial+(i*(A.final_val-A.initial))/A.N), params);
		
		integral_val += (stepsize/6.0) *
		((*pleft_step) + (*pright_step) + 
		(4.0*(*func)(((A.initial+(i*(A.final_val-A.initial))/A.N) +
		(A.initial + (((i-1)*(A.final_val-A.initial))/A.N)))/2.0,params)));
		
		(*pleft_step) = (*pright_step);
		
		//printf("\ndumme summe = %6.10lf\n", h);
	}
	printf("simpson integral = %+6.10lf\n\n", integral_val);
}

