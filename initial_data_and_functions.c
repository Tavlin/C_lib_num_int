#include <math.h>
#include "initial_data_and_functions.h"

FunctionParams g_p_init(double x, double y)
/* initialise the gaussian parameters 
*/
{
	FunctionParams z;
	z.mu = x;
	z.sigma = y;
	return z;
}

InitialData initialdata_init(double a, double b, double c)
{
	InitialData x;
	x.N = a;
	x.initial = b;
	x.final_val = c;
	return x;
}

double gaussian(double x, FunctionParams params)
{
	return (1/(sqrt(2*M_PI)*params.sigma))*
	exp(-pow(x-params.mu, 2)/(2*pow(params.sigma, 2)));
}

double strange_cos(double x, FunctionParams params)
{
	return(x*pow(cos(2*M_PI*pow(x,2)),2));
}
