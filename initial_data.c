#include"initial_data.h"

/*
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

*/
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
