#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "integrals.h"


int main (void)
{

	char* function_name_list[4] = 
	{"Gaussian", "x*cos²(2*pi*x²)", "exp^(-x²)", "1/sqrt(x)"};
	
	char* integral_name_list[7] = 
	{"Left Riemann Sum", "Right Riemann Sum", "Trapezodial Rule",
	"Simpson's Rule", "Trapezodial Rule with semi adaptive stepsize",
	"Midpoint Rule with semiadatipve stepsize",
	"Midpoint Rule with semiadatipve stepsize and open bondary"};
	
	double (*function_list[4])(double x, FunctionParams params) =
	{gaussian, strange_cos, exp_minus_x_sq, rev_sqrt};
	
	double (* integral_list[7])(InitialData A, FunctionParams params,
	double(*func)(double, FunctionParams), double eps) =
	{left_riemann_sum, right_riemann_sum, trapezodial_integral, simpson_integral,
	trapezodial_integral_sas, midpoint_int, midpoint_int_to_inf};
	
	/*
	function_name_list[0] = "Gaussian";
	function_name_list[1] = "x*cos²(2*pi*x²)";
	function_name_list[2] = "exp^(-x²)";
	function_name_list[3] = "1/sqrt(x)";
	*/
	
	
	int check;
	printf("Enter mode (1 for [a,b] or 2 for [a,inf)): ");
		scanf("%i", &check);
		
	if(check == 1)
	{
		double N,a,b,c,mu,sigma;
		printf("Enter stepsize: ");
		scanf("%lf", &N);

		printf("\nEnter lower bound: ");
		scanf("%lf", &a);
		
		printf("\nEnter upper bound: ");
		scanf("%lf", &b);
		
		printf("\nEnter mu: ");
		scanf("%lf", &mu);

		printf("\nEnter sigma: ");
		scanf("%lf", &sigma);

		InitialData A;
		InitialData* pA;
		pA = &A;
		*pA = initialdata_init(N,a,b,c);
	
		FunctionParams norm;
		FunctionParams * pnorm;
		pnorm = &norm;
		*pnorm = g_p_init(mu,sigma); 

		
		printer("Gaussian", "Left Riemann-Sum",
		left_riemann_sum(*pA, *pnorm, gaussian, 0));
		
		printer("Gaussian", "Right Riemann-Sum",
		right_riemann_sum(*pA, *pnorm, gaussian, 0));
		
		printer("Gaussian", "Trapezodial Integral",
		trapezodial_integral(*pA, *pnorm, gaussian, 0));

		printer("Gaussian", "Simpson Integral",
		simpson_integral(*pA, *pnorm, gaussian, 0));
		
		printer("x*cos²(2*pi*x²)", "Left Riemann-Sum",
		left_riemann_sum(*pA, *pnorm, strange_cos, 0));
		
		printer("x*cos²(2*pi*x²)", "Right Riemann-Sum",
		right_riemann_sum(*pA, *pnorm, strange_cos, 0));
		
		printer("x*cos²(2*pi*x²)", "Trapezodial Integral",
		trapezodial_integral(*pA, *pnorm, strange_cos, 0));

		printer("x*cos²(2*pi*x²)", "Simpson Integral",
		simpson_integral(*pA, *pnorm, strange_cos, 0));

		printer("Gaussian", "semiadatipve Trapezodial Integral",
		trapezodial_integral_sas(*pA, *pnorm, gaussian, 0.0000001));
		
		printer("x*cos²(2*pi*x²)", "semiadatipve Trapezodial Integral",
		trapezodial_integral_sas(*pA, *pnorm, strange_cos, 0.0000001));
	
		printer("Gaussian", "semiadatipve Trapezodial Integral",
		midpoint_int(*pA, *pnorm, gaussian, 0.0000001));
		
		printer("x*cos²(2*pi*x²)", "semiadatipve Trapezodial Integral",
		midpoint_int(*pA, *pnorm, strange_cos, 0.0000001));
		
		printer("exp^(-x²)", "Midpointrule with open boundary to infinite", 
		midpoint_int_to_inf(*pA, *pnorm, exp_minus_x_sq, 0.0000001));
		
	}
	
	else if(check == 2)
	{
		double N,a,b,c,mu,sigma;
		printf("Enter stepsize: ");
		scanf("%lf", &N);

		printf("\nEnter lower bound: ");
		scanf("%lf", &a);
		
		printf("\nEnter mu: ");
		scanf("%lf", &mu);

		printf("\nEnter sigma: ");
		scanf("%lf", &sigma);

		InitialData A;
		InitialData* pA;
		pA = &A;
		*pA = initialdata_init(N,a,b,c);
	
		FunctionParams norm;
		FunctionParams * pnorm;
		pnorm = &norm;
		*pnorm = g_p_init(mu,sigma);
		
		
		printer("exp^(-x²)", "Midpointrule with open boundary to infinite", 
		midpoint_int_to_inf(*pA, *pnorm, exp_minus_x_sq, 0.00001));
	}
	return 0;
}
