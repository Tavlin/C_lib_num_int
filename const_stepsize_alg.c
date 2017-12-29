#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "initial_data.h"
#include "integrals.h"

double gaussian(double x, Gaussian params)
{
	return (1/(sqrt(2*M_PI)*params.sigma))*
	exp(-pow(x-params.mu, 2)/(2*pow(params.sigma, 2)));
}


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
	
	Gaussian norm;
	Gaussian * pnorm;
	pnorm = &norm;
	*pnorm = g_p_init(mu,sigma); 

	
	left_riemann_sum(*pA, *pnorm, gaussian);
	right_riemann_sum(*pA, *pnorm, gaussian);
	return 0;
}
