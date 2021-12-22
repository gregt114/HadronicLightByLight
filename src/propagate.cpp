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

    const int samples = int(5e6);     // Number of samples used in integration
    const double cutoff = 0.32;       // Upper bound on Q1 and Q2    


    // Percent uncertainties in parameters a,b,c,d,e,gamma
    const double percent1 = 0.01;
    const double percent2 = 0.005;


    // Parameters for form factor (vary by 0.5% and 0.25%)
    struct params p0 = {a0, b0, c0, d0, e0, gamma_pi0};

    // + percent1
    struct params pa1 = {a0*(1+percent1), b0, c0, d0, e0, gamma_pi0};
    struct params pb1 = {a0, b0*(1+percent1), c0, d0, e0, gamma_pi0};
    struct params pc1 = {a0, b0, c0*(1+percent1), d0, e0, gamma_pi0};
    struct params pd1 = {a0, b0, c0, d0*(1+percent1), e0, gamma_pi0};
    struct params pe1 = {a0, b0, c0, d0, e0*(1+percent1), gamma_pi0};
    struct params pg1 = {a0, b0, c0, d0, e0, gamma_pi0*(1+percent1)};

    // - percent1
    struct params pa2 = {a0*(1-percent1), b0, c0, d0, e0, gamma_pi0};
    struct params pb2 = {a0, b0*(1-percent1), c0, d0, e0, gamma_pi0};
    struct params pc2 = {a0, b0, c0*(1-percent1), d0, e0, gamma_pi0};
    struct params pd2 = {a0, b0, c0, d0*(1-percent1), e0, gamma_pi0};
    struct params pe2 = {a0, b0, c0, d0, e0*(1-percent1), gamma_pi0};
    struct params pg2 = {a0, b0, c0, d0, e0, gamma_pi0*(1-percent1)};

    // + percent2
    struct params pa3 = {a0*(1+percent2), b0, c0, d0, e0, gamma_pi0};
    struct params pb3 = {a0, b0*(1+percent2), c0, d0, e0, gamma_pi0};
    struct params pc3 = {a0, b0, c0*(1+percent2), d0, e0, gamma_pi0};
    struct params pd3 = {a0, b0, c0, d0*(1+percent2), e0, gamma_pi0};
    struct params pe3 = {a0, b0, c0, d0, e0*(1+percent2), gamma_pi0};
    struct params pg3 = {a0, b0, c0, d0, e0, gamma_pi0*(1+percent2)};

    // - percent2
    struct params pa4 = {a0*(1-percent2), b0, c0, d0, e0, gamma_pi0};
    struct params pb4 = {a0, b0*(1-percent2), c0, d0, e0, gamma_pi0};
    struct params pc4 = {a0, b0, c0*(1-percent2), d0, e0, gamma_pi0};
    struct params pd4 = {a0, b0, c0, d0*(1-percent2), e0, gamma_pi0};
    struct params pe4 = {a0, b0, c0, d0, e0*(1-percent2), gamma_pi0};
    struct params pg4 = {a0, b0, c0, d0, e0, gamma_pi0*(1-percent2)};
    

    
    std::cout << "Calculating many integrals (may take a while...) " << std::endl << std::endl;

    // Values of integrals varying given parameter by given percent (takes ~ 5 min)
    double val_a1 = mcIntegrate(&pa1, samples, cutoff);
    double val_b1 = mcIntegrate(&pb1, samples, cutoff);
    double val_c1 = mcIntegrate(&pc1, samples, cutoff);
    double val_d1 = mcIntegrate(&pd1, samples, cutoff);
    double val_e1 = mcIntegrate(&pe1, samples, cutoff);
    double val_g1 = mcIntegrate(&pg1, samples, cutoff);

    double val_a2 = mcIntegrate(&pa2, samples, cutoff);
    double val_b2 = mcIntegrate(&pb2, samples, cutoff);
    double val_c2 = mcIntegrate(&pc2, samples, cutoff);
    double val_d2 = mcIntegrate(&pd2, samples, cutoff);
    double val_e2 = mcIntegrate(&pe2, samples, cutoff);
    double val_g2 = mcIntegrate(&pg2, samples, cutoff);

    double val_a3 = mcIntegrate(&pa3, samples, cutoff);
    double val_b3 = mcIntegrate(&pb3, samples, cutoff);
    double val_c3 = mcIntegrate(&pc3, samples, cutoff);
    double val_d3 = mcIntegrate(&pd3, samples, cutoff);
    double val_e3 = mcIntegrate(&pe3, samples, cutoff);
    double val_g3 = mcIntegrate(&pg3, samples, cutoff);

    double val_a4 = mcIntegrate(&pa4, samples, cutoff);
    double val_b4 = mcIntegrate(&pb4, samples, cutoff);
    double val_c4 = mcIntegrate(&pc4, samples, cutoff);
    double val_d4 = mcIntegrate(&pd4, samples, cutoff);
    double val_e4 = mcIntegrate(&pe4, samples, cutoff);
    double val_g4 = mcIntegrate(&pg4, samples, cutoff);

    std::cout << "Integration up to Q < " << cutoff << std::endl << std::endl;

    std::cout.precision(17);    // Decimal precision for printed output

    // Partial derivatives using double sided finite difference
    std::cout << "Partials - parameters varied by " << 100*percent1 << "\%\n";
    std::cout << "a     :       " << (val_a1 - val_a2)/(2*percent1*a0) << std::endl;
    std::cout << "b     :       " << (val_b1 - val_b2)/(2*percent1*b0) << std::endl;
    std::cout << "c     :       " << (val_c1 - val_c2)/(2*percent1*c0) << std::endl;
    std::cout << "d     :       " << (val_d1 - val_d2)/(2*percent1*d0) << std::endl;
    std::cout << "e     :       " << (val_e1 - val_e2)/(2*percent1*e0) << std::endl;
    std::cout << "gamma :       " << (val_g1 - val_g2)/(2*percent1*gamma_pi0) << std::endl << std::endl;

    std::cout << "Partials - parameters varied by " << 100*percent2 << "\%\n";
    std::cout << "a     :       " << (val_a3 - val_a4)/(2*percent2*a0) << std::endl;
    std::cout << "b     :       " << (val_b3 - val_b4)/(2*percent2*b0) << std::endl;
    std::cout << "c     :       " << (val_c3 - val_c4)/(2*percent2*c0) << std::endl;
    std::cout << "d     :       " << (val_d3 - val_d4)/(2*percent2*d0) << std::endl;
    std::cout << "e     :       " << (val_e3 - val_e4)/(2*percent2*e0) << std::endl;
    std::cout << "gamma :       " << (val_g3 - val_g4)/(2*percent2*gamma_pi0) << std::endl << std::endl;



    std::cout << std::endl;

    return 0;
}


// Calculating many integrals (may take a while...) 

// Integration up to Q < 0.32

// Partials - parameters varied by 2%
// a     :       -3.3255397848918805e-11
// b     :       2.5390081964826805e-12
// c     :       9.6649245908106711e-13
// d     :       4.5293740842159872e-13
// e     :       1.1512917587814967e-13
// gamma :       0.028049832503515056

// Partials - parameters varied by 1%
// a     :       -3.3720972009572263e-11
// b     :       2.5927594821263501e-12
// c     :       9.8838045327526882e-13
// d     :       4.8949333115808346e-13
// e     :       1.7264069732678279e-13
// gamma :       0.02804983250351514


