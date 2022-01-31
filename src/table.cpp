#include <iostream>
#include <cmath>
#include <chrono>


#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

#include "../include/functions.h"


// First function to integrate (integrand1)
double f1(double x[], size_t dim, void * p)
{
    (void)(dim);
    (void)(p);

    // By default uses LMD+V form factor. Change this to:
    // return integrand1_q6(x[0], x[1], x[2])
    // if you want to use the Q^6 expansion instead
    return integrand1(x[0], x[1], x[2]);
}

// Second function to integrate (integrand2)
double f2(double x[], size_t dim, void * p)
{
    (void)(dim);
    (void)(p);

    // By default uses LMD+V form factor. Change this to:
    // return integrand2_q6(x[0], x[1], x[2])
    // if you want to use the Q^6 expansion instead
    return integrand2(x[0], x[1], x[2]);
}




int main()
{ 
  const double limits[] = {0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 5.0, 10.0, 20.0};                // Momentum cutoffs for integration on Q1 and Q2

  double res1, err1, res2, err2 = 0;

  std::cout << "Cutoff\t\tValue\n\n";

  // Iterate through each cutoff and calculate a_mu
  for(int i=0; i < sizeof(limits)/sizeof(limits[0]); i++){
      
    double limit = limits[i];

    // Number of samples used in Monte Carlo integration
    size_t calls = int(10e6);
    
    // Limits of integration
    double xl[3] = { 0, 0, -1 };
    double xu[3] = { limit, limit, 1 };

    // MISER Setup
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_monte_function G1 = { &f1, 3, 0 };
    gsl_monte_function G2 = { &f2, 3, 0 };

    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    // Perform calculations
    {
    gsl_monte_miser_state *s = gsl_monte_miser_alloc (3);
    gsl_monte_miser_integrate (&G1, xl, xu, 3, calls, r, s, &res1, &err1);
    gsl_monte_miser_integrate (&G2, xl, xu, 3, calls, r, s, &res2, &err2);
    gsl_monte_miser_free (s);
    }
    gsl_rng_free (r);

    // Print result
    std::cout << limit << "\t\t" << pow(alpha/M_PI, 3) * (res1 + res2) << "\n";
  }   


  return 0;
}


