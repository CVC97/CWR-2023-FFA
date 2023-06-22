#ifndef RNG_H
#define RNG_H

#include "cvc_numerics.h"


// Tuple of 2 gaussian-distributed random numbers with static random-number-generator: Polarmethode
struct cvc_tuple_2 cvc_random_gaussian(gsl_rng* generator);

// Volumen eines Hyperquaders R
double domain_volume(int D, double R[]);

// 2-Dimensional Integration: Monte-Carlo
double cvc_mc_integrate_2D(gsl_rng* generator, int A(double, double), double a_x, double b_x, double a_y, double b_y, int N, double f(double, double));

// MC-Berechnung der Dichte im Hyperquader R [xi_min, xi_max] mit D Dimensionen für die Dichtefunktion integrand()  
double cvc_mc_integrate(gsl_rng* generator, int D, double integrand(int, double*), double R[], int N);

#endif