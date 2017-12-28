#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

struct gaussian_parameters
{
	double mu;
	double sigma;
};

struct gaussian_parameters g_p_init(double x, double y)
/* initialise the gaussian parameters 
*/
{

	struct gaussian_parameters v;
	v.mu = x;
	v.sigma = y;
	return v;
}


double gaussian(double x, struct gaussian_parameters params)
{
	params.mu = 0;
	params.sigma = 1;
	return (1/(sqrt(2*M_PI)*params.sigma))*
	exp(-pow(x-params.mu, 2)/(2*pow(params.sigma, 2)));
}

//left Riemann sum
double left_riemann_sum(double N, double initial, double final, struct gaussian_parameters params,
double(*func)(double, struct gaussian_parameters)) 
/*initial instead of init, cuz __init__ is SPECIAL, at least somewhere else :)
*/
{
	double left_sum = 0;
	for(double i = 1; i <= N; i++)
	{
		left_sum += ((initial+(i*(final-initial))/N) - /*upper bound*/
		(initial + (((i-1)*(final-initial))/N))) * /*lower bound*/
		(*func)((initial + (((i-1)*(final-initial))/N)), params);
	}
	printf("left riemann sum = %lf\n", left_sum);
}


//right Rieman sum
double right_riemann_sum(double N, double initial, double final, struct gaussian_parameters params,
double(*func)(double, struct gaussian_parameters))
{
	double right_sum = 0;
	for(double i = 1; i <= N; i++)
	{
		right_sum += ((initial+(i*(final-initial))/N) - /*upper bound*/
		(initial + (((i-1)*(final-initial))/N))) * /*lower bound*/
		(*func)(((initial+(i*(final-initial))/N)), params);
	}
	printf("right riemann sum = %lf\n", right_sum);
}






int main (void)
{
	double N,a,b;
	printf("Enter stepsize: ");
	scanf("%lf", &N);
	a = -1;
	b = 1;
	struct gaussian_parameters normaldistribution;
	normaldistribution = g_p_init(0,1); 
	left_riemann_sum(N, a, b, normaldistribution, gaussian);
	right_riemann_sum(N, a, b, normaldistribution, gaussian);
	return 0;
}
