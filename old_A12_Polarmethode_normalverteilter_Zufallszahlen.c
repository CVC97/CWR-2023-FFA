#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_rng.h>


// Struktur eines 2er Tupels
struct tuple_2 {
    double x1;
    double x2;
};


// Returniert 2 normalverteile Zuvallsgrößen für gegebenen (GSL_RNG) Zufallsgenerator
struct tuple_2 random_gaussian(gsl_rng* generator) {
    // Seed with time applied to given generator
    gsl_rng_set(generator, time(NULL));     
    for (int i = 0; i < 10; i++) {
        printf("%.4f ", gsl_rng_uniform(generator));
    }
    printf("\n");

    // Erstelung 2 Zufallszahlen u, v
    double u, v, r = 0, m;
    while (r > 1 || r == 0){
        u = (gsl_rng_uniform(generator) * 2) - 1;
        v = (gsl_rng_uniform(generator) * 2) - 1;
        r = pow(u, 2) + pow(v, 2);
    }
    m = sqrt(-2*log(r)/r);

    // Tupel mit 2 normalverteilten Zufallszahlen
    struct tuple_2 random_2;
    random_2.x1 = u*m;
    random_2.x2 = v*m;
    return random_2;
}


int main(void) {
    // GSL's Taus generator:
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);

    random_gaussian(rng);
    gsl_rng_free(rng);
    return 0;
}