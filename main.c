#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "integrals.h"


int main (void)
{
	int check;
	printf("1 for [a,b]\n2 for [a,inf):");
		scanf("%i", &check);
		
	if(check == 1)
	{
		double N,a,b,mu,sigma;
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
		*pA = initialdata_init(N,a,b);
	
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
		midpoint_int_to_inf(*pA, *pnorm, exp_minus_x_sq, 0.00001));
		
	}
	
	else if(check == 2)
	{
		double N,a,b,mu,sigma;
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
		*pA = initialdata_init(N,a,b);
	
		FunctionParams norm;
		FunctionParams * pnorm;
		pnorm = &norm;
		*pnorm = g_p_init(mu,sigma);
		
		
		printer("exp^(-x²)", "Midpointrule with open boundary to infinite", 
		midpoint_int_to_inf(*pA, *pnorm, exp_minus_x_sq, 0.00001));
	}
	return 0;
}
