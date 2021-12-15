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
  double const limit = 20;                // Momentum cutoff for integration on Q1 and Q2


  double res, err, res2, err2;

  // Limits of integration
  double xl[3] = { 0, 0, -1 };
  double xu[3] = { limit, limit, 1 };

  // Number of samples used in Monte Carlo integration
  size_t calls = int(4e7);


  // Setup for GSL Monte Carlo algorithm
  const gsl_rng_type *T;
  gsl_rng *r;
  gsl_monte_function G = { &f1, 3, 0 };
  gsl_monte_function G2 = { &f2, 3, 0 };
  gsl_rng_env_setup ();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);


  // VEGAS Algorithm
  {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);

    gsl_monte_vegas_integrate (&G, xl, xu, 3, 10000, r, s,
                               &res, &err);

    gsl_monte_vegas_integrate (&G2, xl, xu, 3, 10000, r, s,
                            &res2, &err2);

    printf ("converging...\n");

    do
      {
        gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/5, r, s,
                                   &res, &err);
        gsl_monte_vegas_integrate (&G2, xl, xu, 3, calls/5, r, s,
                                   &res2, &err2);
      }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);


  // Print results
  std::cout.precision(17);
  std::cout << "Integral 1:\t" << res << "\tSigma: " <<  err << "\n";
  std::cout << "Integral 2:\t" << res2 << "\tSigma: " <<  err2 << "\n";
  std::cout << "Final Result:\t" << pow(alpha/M_PI, 3) * (res + res2) << "\n\n";

  gsl_monte_vegas_free (s);
  }


  gsl_rng_free (r);

  return 0;
}


/*
0.0489815	0.000107344
0.00130036	3.6389e-06
converging...
Integral 1:	0.048916	Sigma: 6.3488e-08
Integral 2:	0.00130896	Sigma: 2.22985e-09
Final Result:	6.29457e-10

400 mil samples
VEGAS Algo
*/


