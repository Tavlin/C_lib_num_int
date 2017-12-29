struct gaussian_parameters
{
	double mu;
	double sigma;
};

struct bounds_stepsize
{
	double N;
	double initial;
	double final;
};


// structure for gaussian parameters
typedef struct gaussian_parameters Gaussian;

//structure for the integral-bounds and the stepsize
typedef struct bounds_stepsize InitialData;


//declaration
Gaussian g_p_init(double x, double y);

//declaration
InitialData initialdata_init(double a, double b, double c);
