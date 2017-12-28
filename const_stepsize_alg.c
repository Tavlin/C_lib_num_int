#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

struct gaussian_parameters
{
	double mu;
	double sigma;
};

double gaussian(double x, struct gaussian_parameters params)
{
	params.mu = 0;
	params.sigma = 1;
	return (1/(sqrt(2*M_PI)*params.sigma))*
	exp(-pow(x-params.mu, 2)/(2*pow(params.sigma, 2)));
}

//left Riemann sum
double left_riemann_sum(double N, double initial, double final, void* params,
double(*func)(double, void*), double h) 
/*initial instead of init, cuz __init__ is SPECIAL, at least somewhere else :)
*/
{
	double left_sum = 0;
	for(double i = 1; i <= N; i++)
	{
		left_sum += ((initial+(i*(final-initial))/N) - /*upper bound*/
		(initial + (((i-1)*(final-initial))/N))) * /*lower bound*/
		(*func)((initial + (((i-1)*(final-initial))/N)),params);
	}
}


//right Rieman sum
double right_riemann_sum(double N, double initial, double final, void* params,
double(*func)(double, void*))
{
	double right_sum = 0;
	for(double i = 1; i <= N; i++)
	{
		right_sum += ((initial+(i*(final-initial))/N) - /*upper bound*/
		(initial + (((i-1)*(final-initial))/N))) * /*lower bound*/
		(*func)(((initial+(i*(final-initial))/N)),params);
	}
}






int main (void)
{
	return 0;
}
