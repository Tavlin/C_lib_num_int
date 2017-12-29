#include <math.h>
#include "initial_data_and_functions.h"

Gaussian g_p_init(double x, double y)
/* initialise the gaussian parameters 
*/
{
	Gaussian z;
	z.mu = x;
	z.sigma = y;
	return z;
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
