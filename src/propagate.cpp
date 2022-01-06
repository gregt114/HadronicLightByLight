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
    double gamma;
};

// Function to integrate (integrand1)
double func1 (double x[], size_t dim, void * p) {
   struct params * fp = (struct params *)p;
   return integrand1_q6(x[0], x[1], x[2], fp->a, fp->b, fp->c, fp->d, fp->e, fp->gamma);
}

// Function to integrate (integrand2)
double func2 (double x[], size_t dim, void * p) {
   struct params * fp = (struct params *)p;
   return integrand2_q6(x[0], x[1], x[2], fp->a, fp->b, fp->c, fp->d, fp->e, fp->gamma);
}



/*
Function that calculates a^(HLbL) using Monte Carlo integration and the given parameters (uses Q^6 approximation)

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

    const int samples = int(10e6);     // Number of samples used in integration
    const double cutoff = 0.32;        // Upper bound on Q1 and Q2    


    std::cout << "Calculating many integrals (may take a while...) " << std::endl << std::endl;
    std::cout.precision(15);

    const double percents[] = {0.01, 0.02};             // Percent uncertainties in parameters a,b,c,d,e,gamma
    char paramNames[] = {'a', 'b', 'c', 'd', 'e', 'g'}; // Labels for parameters
    double avgs[] = {a0, b0, c0, d0, e0, gamma_pi0};    // Mean values for parameters

    for(unsigned i=0; i < sizeof(percents) / sizeof(percents[0]); i++){ // Iterate over each % in percents[]
        
        double percent = percents[i];
        std::cout << "Partials - parameters varied by " << 100*percent << "\%\n";

        struct params plus_p[] = {  // + percent
            {a0*(1+percent), b0, c0, d0, e0, gamma_pi0},
            {a0, b0*(1+percent), c0, d0, e0, gamma_pi0},
            {a0, b0, c0*(1+percent), d0, e0, gamma_pi0},
            {a0, b0, c0, d0*(1+percent), e0, gamma_pi0},
            {a0, b0, c0, d0, e0*(1+percent), gamma_pi0},
            {a0, b0, c0, d0, e0, gamma_pi0*(1+percent)}
        };

        struct params minus_p[] = {  // - percent
            {a0*(1-percent), b0, c0, d0, e0, gamma_pi0},
            {a0, b0*(1-percent), c0, d0, e0, gamma_pi0},
            {a0, b0, c0*(1-percent), d0, e0, gamma_pi0},
            {a0, b0, c0, d0*(1-percent), e0, gamma_pi0},
            {a0, b0, c0, d0, e0*(1-percent), gamma_pi0},
            {a0, b0, c0, d0, e0, gamma_pi0*(1-percent)}
        };

        for(unsigned j=0; j < sizeof(plus_p) / sizeof(plus_p[0]); j++){ // Iterate over each parameter for given %
            double val1 = mcIntegrate(&plus_p[j], samples, cutoff);
            double val2 = mcIntegrate(&minus_p[j], samples, cutoff);
            std::cout << paramNames[j]  << "\t:\t"<< (val1 - val2)/(2*percent*avgs[j]) << std::endl;
        }
        std::cout << std::endl << std::endl;
    }
    

    return 0;
}
