#ifndef RNG_H
#define RNG_H

#include "cvc_numerics.h"


// Tuple of 2 gaussian-distributed random numbers with static random-number-generator: Polarmethode
struct cvc_tuple_2 cvc_random_gaussian(gsl_rng* generator);

// MC-Berechnung der Dichte im Hyperquader [ai, bi] mit dim Dimensionen für die Dichtefunktion func()  
double cvc_mc_integrate(double func(double*, int), double a[], double b[], int dim, int N, gsl_rng* generator);

// 2-Dimensional Integration: Monte-Carlo
double cvc_mc_integrate_2D(int A(double, double), double a_x, double b_x, double a_y, double b_y, int N, double f(double, double), gsl_rng* generator);

#endif