#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "integrals.h"

int main (void)
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

	left_riemann_sum(*pA, *pnorm, gaussian);
	right_riemann_sum(*pA, *pnorm, gaussian);
	trapezodial_integral(*pA, *pnorm, gaussian);
	simpson_integral(*pA, *pnorm, gaussian);
	
	left_riemann_sum(*pA, *pnorm, strange_cos);
	right_riemann_sum(*pA, *pnorm, strange_cos);
	trapezodial_integral(*pA, *pnorm, strange_cos);
	simpson_integral(*pA, *pnorm, strange_cos);
	
	return 0;
}
