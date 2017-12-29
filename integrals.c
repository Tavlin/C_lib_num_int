#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "initial_data.h"
#include "integrals.h"

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
	printf("left riemann sum = %+6.10lf\n", left_sum);
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
	printf("right riemann sum = %+6.10lf\n", right_sum);
}