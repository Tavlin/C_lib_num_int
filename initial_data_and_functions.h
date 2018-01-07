// Definition of function parameters
struct function_parameters
{
	double mu;
	double sigma;
};

// Definition for the integral-bounds and the stepsize
struct bounds_stepsize
{
	double N;
	double initial;
	double final_val;
};


// structure for function parameters
typedef struct function_parameters FunctionParams;

// structure for the integral-bounds and the stepsize
typedef struct bounds_stepsize InitialData;


// declaration
FunctionParams g_p_init(double x, double y);

//d eclaration
InitialData initialdata_init(double a, double b, double c);

double gaussian(double x, FunctionParams params);

double strange_cos(double x, FunctionParams params);
