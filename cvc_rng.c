#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_rng.h>
#include "cvc_numerics.h"


// Tupel aus 2 normalverteile Zuvallsgrößen für gegebenen (GSL_RNG) Zufallsgenerator: Polarmethode
struct cvc_tuple_2 cvc_random_gaussian(gsl_rng* generator) {   
    double u, v, r = 0, m;                                          // Erstellung 2 Zufallszahlen u, v
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


// Volumen eines Hyperquaders R [xi_min, xi_max] mit D Dimensionen
double domain_volume(int D, double R[]) {
    double volume_R = 1;
    for (int i_dim = 0; i_dim < D; i_dim++) {
        volume_R *= (R[2*i_dim + 1] - R[2*i_dim]);
    }
    return volume_R;
}


// 2-Dimensionale Integration: MC
double cvc_mc_integrate_2D(gsl_rng* generator, int A(double, double), double a_x, double b_x, double a_y, double b_y, int N, double f(double, double)) {
    double x, y, sum = 0, R_area = (b_x - a_x) * (b_y - a_y);
    for (int i = 0; i < N; i++) {
        x = gsl_rng_uniform(generator) * (b_x - a_x) + a_x;
        y = gsl_rng_uniform(generator) * (b_y - a_y) + a_y;
        sum += f(x, y) * A(x, y);
    }
    return R_area * sum / N;
}


// MC-Berechnung der Dichte im Hyperquader R [xi_min, xi_max] mit D Dimensionen für die Dichtefunktion integrand()  
double cvc_mc_integrate(gsl_rng* generator, int D, double integrand(int, double*), double R[], int N) { 
    double total_volume = domain_volume(D, R);                      // Gesamtvolumen
    double density_sum = 0;                                         // Laufvariable der Dichtebestimmung
    double *x = (double*) malloc(sizeof(double) * D);               // Erstellung N Zufallsvektoren x[]
    for (int i = 0; i < N; i++) {
        for (int i_dim = 0; i_dim < D; i_dim++) {
            x[i_dim] = gsl_rng_uniform(generator) * (R[2*i_dim + 1] - R[2*i_dim]) + R[2*i_dim];
        }
        density_sum += integrand(D, x);
    }
    free(x);
    return total_volume / N * density_sum;
}