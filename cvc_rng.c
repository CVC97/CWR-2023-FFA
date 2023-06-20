#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_rng.h>
#include "cvc_numerics.h"


// Tupel aus 2 normalverteile Zuvallsgrößen für gegebenen (GSL_RNG) Zufallsgenerator: Polarmethode
struct cvc_tuple_2 cvc_random_gaussian(gsl_rng* generator) {
    gsl_rng_set(generator, time(NULL));    
    double u, v, r = 0, m;                                          // Erstelung 2 Zufallszahlen u, v
    while (r > 1 || r == 0){
        u = (gsl_rng_uniform(generator) * 2) - 1;
        v = (gsl_rng_uniform(generator) * 2) - 1;
        r = pow(u, 2) + pow(v, 2);
    }
    m = sqrt(-2*log(r)/r);
    struct cvc_tuple_2 random_2;                                    // Tupel mit 2 normalverteilten Zufallszahlen
    random_2.x1 = u*m;
    random_2.x2 = v*m;
    return random_2;
}  


// MC-Berechnung der Dichte im Hyperquader [ai, bi] mit dim Dimensionen für die Dichtefunktion func()  
double cvc_mc_integrate(gsl_rng* generator, double func(double*, int), double a[], double b[], int dim, int N) {
    gsl_rng_set(generator, time(NULL));  
    double V = 1, density_sum = 0;                                  // Gesamtvolumen
    for (int i_dim = 0; i_dim < dim; i_dim++) {
        V *= b[i_dim] - a[i_dim];
    }
    double *y = (double*) malloc(sizeof(double) * dim);             // Erstellung N Zufallsvektoren y[]
    for (int i = 0; i < N; i++) {
        for (int i_dim = 0; i_dim < dim; i_dim++) {
            *(y + i_dim) = gsl_rng_uniform(generator) * (b[i_dim] - a[i_dim]) + a[i_dim];
        }
        density_sum += func(y, dim);
    }
    free(y);
    return V / N * density_sum;
}


// 2-Dimensionale Integration: MC
double cvc_mc_integrate_2D(gsl_rng* generator, int A(double, double), double a_x, double b_x, double a_y, double b_y, int N, double f(double, double)) {
    gsl_rng_set(generator, time(NULL)); 
    double x, y, sum = 0, R_area = (b_x - a_x) * (b_y - a_y);
    for (int i = 0; i < N; i++) {
        x = gsl_rng_uniform(generator) * (b_x - a_x) + a_x;
        y = gsl_rng_uniform(generator) * (b_y - a_y) + a_y;
        sum += f(x, y) * A(x, y);
    }
    return R_area * sum / N;
}