#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct gaussian_parameters
{
	double mu;
	double sigma;
}Gaussian;

typedef struct bounds_stepsize
{
	double N;
	double initial;
	double final;
}InitialData;

Gaussian g_p_init(double x, double y)
/* initialise the gaussian parameters 
*/
{

	Gaussian v;
	v.mu = x;
	v.sigma = y;
	return v;
}

InitialData initialdata_init(double a, double b, double c)
{
	InitialData x;
	x.N = a;
	x.initial = b;
	x.final = c;
	return x;
}

double gaussian(double x, Gaussian params)
{
	return (1/(sqrt(2*M_PI)*params.sigma))*
	exp(-pow(x-params.mu, 2)/(2*pow(params.sigma, 2)));
}

//left Riemann sum
double left_riemann_sum(InitialData A, Gaussian params,
double(*func)(double, Gaussian)) 
/*initial instead of init, cuz __init__ is SPECIAL, at least somewhere else :)
*/
{
	double left_sum = 0;
	for(double i = 1; i <= A.N; i++)
	{
		left_sum += ((A.initial+(i*(A.final-A.initial))/A.N) - /*upper bound*/
		(A.initial + (((i-1)*(A.final-A.initial))/A.N))) * /*lower bound*/
		(*func)((A.initial + (((i-1)*(A.final-A.initial))/A.N)), params);
	}
	printf("left riemann sum = %lf\n", left_sum);
}


//right Rieman sum
double right_riemann_sum(InitialData A, Gaussian params,
double(*func)(double, Gaussian)) 
{
	double right_sum = 0;
	for(double i = 1; i <= A.N; i++)
	{
		right_sum += ((A.initial+(i*(A.final-A.initial))/A.N) - /*upper bound*/
		(A.initial + (((i-1)*(A.final-A.initial))/A.N))) * /*lower bound*/
		(*func)(((A.initial+(i*(A.final-A.initial))/A.N)), params);
	}
	printf("right riemann sum = %lf\n", right_sum);
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
	A = initialdata_init(N,a,b);
	
	Gaussian norm;
	norm = g_p_init(mu,sigma); 
	
	left_riemann_sum(A, norm, gaussian);
	right_riemann_sum(A, norm, gaussian);
	return 0;
}
