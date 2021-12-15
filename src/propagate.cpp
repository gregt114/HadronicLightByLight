#include <iostream>
#include <fstream>
#include <cmath>

#include <random>
#include <chrono>

#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

#include "../include/functions.h"



// Defines parameter structure for integration
struct params 
{
    double a;
    double b;
    double c;
    double d;
    double e;
};

// Function to integrate (integrand1)
double func1 (double x[], size_t dim, void * p) {
   struct params * fp = (struct params *)p;
   return integrand1_q6(x[0], x[1], x[2], fp->a, fp->b, fp->c, fp->d, fp->e);
}

// Function to integrate (integrand2)
double func2 (double x[], size_t dim, void * p) {
   struct params * fp = (struct params *)p;
   return integrand2_q6(x[0], x[1], x[2], fp->a, fp->b, fp->c, fp->d, fp->e);
}



/*
Function that calculates a^(HLbL) using Monte Carlo integration and the given parameters

*p      - pointer to params structure that contains a,b,c,d,e values
samples - number of samples used during integration
bound   - the upper limit of integration on Q1 and Q2
*/
double mcIntegrate(struct params *p, int samples, double bound)
{
    // Variables to store results of integration
    double res1, res2;
    double err1, err2;

    // Lower and upper bounds for integration
    double xl[3] = { 0, 0, -1 };
    double xu[3] = { bound, bound, 1 };

    // Setup GSL Monte Carlo Algorithm
    const gsl_rng_type *T;
    gsl_rng *r;

    // Setup GSL Monte Functions
    gsl_monte_function G1, G2;
    G1.f = &func1;
    G1.dim = 3;
    G1.params = p;
    G2.f = &func2;
    G2.dim = 3;
    G2.params = p;

    // More Setup
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    // Perform integration
    {
    gsl_monte_miser_state *s = gsl_monte_miser_alloc (3);
    gsl_monte_miser_integrate (&G1, xl, xu, 3, samples, r, s,
                                &res1, &err1);
    gsl_monte_miser_free (s);
    }

    {
    gsl_monte_miser_state *s = gsl_monte_miser_alloc (3);
    gsl_monte_miser_integrate (&G2, xl, xu, 3, samples, r, s,
                                &res2, &err2);
    gsl_monte_miser_free (s);
    }

    // Return final result
    return pow(alpha/M_PI, 3) * (res1 + res2);
    
}






int main()
{

    int samples = int(5e6);     // Number of samples used in integration
    double cutoff = 0.32;       // Upper bound on Q1 and Q2
    std::cout.precision(17);    // Decimal precision for printed output


    // Percent uncertainties in parameters a,b,c,d,e
    double percent1 = 0.033;
    double percent2 = 0.05;
    double percent3 = 0.1;


    // Parameters for form factor
    struct params p0 = {a0, b0, c0, d0, e0};

    struct params p3 = {a0*(1+percent1), b0, c0, d0, e0};
    struct params p3_2 = {a0*(1-percent1), b0, c0, d0, e0};

    struct params p5 = {a0*(1+percent2), b0, c0, d0, e0};
    struct params p5_2 = {a0*(1-percent2), b0, c0, d0, e0};

    struct params p10 = {a0*(1+percent3), b0, c0, d0, e0};
    struct params p10_2 = {a0*(1-percent3), b0, c0, d0, e0};
    

    // Value of the integral using the mean value for 'a'
    double val0 = mcIntegrate(&p0, samples, cutoff);
    std::cout << std::endl;
    std::cout << "Result with mean value of a                 : " << val0 << std::endl;
    std::cout << "Calculating more integrals (may take a little...) " << std::endl << std::endl;

    // Values of integrals using varied 'a' parameters
    double val1 = mcIntegrate(&p3, samples, cutoff);
    double val2 = mcIntegrate(&p3_2, samples, cutoff);
    
    double val3 = mcIntegrate(&p5, samples, cutoff);
    double val4 = mcIntegrate(&p5_2, samples, cutoff);

    double val5 = mcIntegrate(&p10, samples, cutoff);
    double val6 = mcIntegrate(&p10_2, samples, cutoff);

    
    std::cout << "a -+ 3.3%                                   : " << val1 << " , " << val2 << std::endl;
    std::cout << "a -+ 5%                                     : " << val3 << " , " << val4 << std::endl;
    std::cout << "a -+ 10%                                    : " << val5 << " , " << val6 << std::endl;
    std::cout << std::endl;
    std::cout << "Percent change in final result (a -+ 3.3%)  : " << abs(val1-val0)/val0 * 100 << " , " << abs(val2-val0)/val0 * 100 << std::endl;
    std::cout << "Percent change in final result (a -+ 5%)    : " << abs(val3-val0)/val0 * 100 << " , " << abs(val4-val0)/val0 * 100 << std::endl;
    std::cout << "Percent change in final result (a -+ 10%)   : " << abs(val5-val0)/val0 * 100 << " , " << abs(val6-val0)/val0 * 100 << std::endl;



    std::cout << std::endl;

    return 0;
}
