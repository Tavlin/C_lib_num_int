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

	
		left_riemann_sum(*pA, *pnorm, gaussian, 0);
		right_riemann_sum(*pA, *pnorm, gaussian, 0);
		trapezodial_integral(*pA, *pnorm, gaussian, 0);
		simpson_integral(*pA, *pnorm, gaussian, 0);
	
		left_riemann_sum(*pA, *pnorm, strange_cos, 0);
		right_riemann_sum(*pA, *pnorm, strange_cos, 0);
		trapezodial_integral(*pA, *pnorm, strange_cos, 0);
		simpson_integral(*pA, *pnorm, strange_cos, 0);
	
		trapezodial_integral_sas(*pA, *pnorm, gaussian, 0.0000001);
		trapezodial_integral_sas(*pA, *pnorm, strange_cos, 0.0000001);
		midpoint_int(*pA, *pnorm, gaussian, 0.00000001);
		midpoint_int(*pA, *pnorm, strange_cos, 0.0000001);
		
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
		
		midpoint_int_to_inf(*pA, *pnorm, gaussian, 0.000001);
	}
	return 0;
}
