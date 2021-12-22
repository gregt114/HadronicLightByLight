#include <iostream>
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



// LMD+V Form Factor Integrands
double f1_LMDV(double x[], size_t dim, void * p)
{
    (void)(dim);
    (void)(p);
    return integrand1(x[0], x[1], x[2]);
}

double f2_LMDV(double x[], size_t dim, void * p)
{
    (void)(dim);
    (void)(p);
    return integrand2(x[0], x[1], x[2]);
}

// Integrands Q^4 Expansion
double f1_q4(double x[], size_t dim, void * p)
{
    (void)(dim);
    (void)(p);
    return integrand1_q4(x[0], x[1], x[2]);
}

double f2_q4(double x[], size_t dim, void * p)
{
    (void)(dim);
    (void)(p);
    return integrand2_q4(x[0], x[1], x[2]);
}

// Integrands Q^6 Expansion
double f1_q6(double x[], size_t dim, void * p)
{
    (void)(dim);
    (void)(p);
    return integrand1_q6(x[0], x[1], x[2]);
}

double f2_q6(double x[], size_t dim, void * p)
{
    (void)(dim);
    (void)(p);
    return integrand2_q6(x[0], x[1], x[2]);
}




int main()
{    
    size_t calls = int(40e6);                  // Number of samples used in Monte Carlo integration
    double limit = 0.1; //0.55                // Upper integration bound on Q1 and Q2


    // Variables to store results
    double res1, res2, res3, res4, res5, res6;
    double err1, err2, err3, err4, err5, err6;

    // Lower and upper bounds for integration
    double xl[3] = { 0, 0, -1 };
    double xu[3] = { limit, limit, 1 };



    // Setup GSL Monte Carlo Algorithm
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_monte_function G1 = { &f1_LMDV, 3, 0 };
    gsl_monte_function G2 = { &f2_LMDV, 3, 0 };
    gsl_monte_function G3 = { &f1_q4, 3, 0 };
    gsl_monte_function G4 = { &f2_q4, 3, 0 };
    gsl_monte_function G5 = { &f1_q6, 3, 0 };
    gsl_monte_function G6 = { &f2_q6, 3, 0 };
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    // Perform calculations
    std::cout << "Calculating many integrals (may take a while...) " << std::endl << std::endl;
    {
    gsl_monte_miser_state *s = gsl_monte_miser_alloc (3);
    gsl_monte_miser_integrate (&G1, xl, xu, 3, calls, r, s,
                                &res1, &err1);

    gsl_monte_miser_integrate (&G2, xl, xu, 3, calls, r, s,
                                &res2, &err2);

    // Uncomment this if you also want the accuracy of the Q^4 approximation
    gsl_monte_miser_integrate (&G3, xl, xu, 3, calls, r, s,
                            &res3, &err3);

    gsl_monte_miser_integrate (&G4, xl, xu, 3, calls, r, s,
                                &res4, &err4);

    gsl_monte_miser_integrate (&G5, xl, xu, 3, calls, r, s,
                            &res5, &err5);

    gsl_monte_miser_integrate (&G6, xl, xu, 3, calls, r, s,
                                &res6, &err6);
    gsl_monte_miser_free (s);

    // Calculate final results for each form factor
    double LMDV = pow(alpha/M_PI, 3) * (res1 + res2);
    double Q4 = pow(alpha/M_PI, 3) * (res3 + res4);
    double Q6 = pow(alpha/M_PI, 3) * (res5 + res6);

    std::cout.precision(10);
    std::cout << "Integration to Q < " << limit << "\n";
    std::cout << "Integral 1 (LMD+V):\t" << res1 << "\tSigma: " <<  err1 << "\n";
    std::cout << "Integral 2 (LMD+V):\t" << res2 << "\tSigma: " <<  err2 << "\n";
    std::cout << std::endl;
    // Uncomment this(and the print statements below) if you also want the accuracy of the Q^4 approximation
    std::cout << "Integral 1 (Q4)   :\t" << res3 << "\tSigma: " <<  err3 << "\n";
    std::cout << "Integral 2 (Q4)   :\t" << res4 << "\tSigma: " <<  err4 << "\n";
    std::cout << std::endl;
    std::cout << "Integral 1 (Q6)   :\t" << res5 << "\tSigma: " <<  err5 << "\n";
    std::cout << "Integral 2 (Q6)   :\t" << res6 << "\tSigma: " <<  err6 << "\n";
    std::cout << std::endl;
    std::cout << "Final LMD+V       :\t" << LMDV << "\n";
    std::cout << "Final Q4          :\t" << Q4 << "\n";
    std::cout << "Final Q6          :\t" << Q6 << "\n\n";
    std::cout << std::endl;
    std::cout << "% Error Q4 = " << 100 * abs(LMDV - Q4) / LMDV << "\n";
    std::cout << "% Error Q6 = " << 100 * abs(LMDV - Q6) / LMDV << "\n";
    }




  gsl_rng_free (r);

  return 0;
}







