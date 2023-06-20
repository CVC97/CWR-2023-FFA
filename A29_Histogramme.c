#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_rng.h>
#include "cvc_numerics.h"
#include "cvc_rng.h"


int main(void) {
    int n_bins = 30;
    double rand1, rand2;

    gsl_histogram *hist_with_gsl = gsl_histogram_alloc(n_bins);
    gsl_histogram_set_ranges_uniform(hist_with_gsl, -5, 5);

    // Wahl des Zufallsgenerators: gsl_rng_mt19937
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

    for (int i = 0; i < 5e6; i++) {
        struct cvc_tuple_2 random_tuple_2 = cvc_random_gaussian(rng);
        rand1 = random_tuple_2.x1;
        rand2 = random_tuple_2.x2;

        printf("iter: %i\trand1: %f, rand2: %f\n", i, rand1, rand2);

        gsl_histogram_increment(hist_with_gsl, rand1);
        gsl_histogram_increment(hist_with_gsl, rand2);
    }

    FILE* file_histogram = fopen("data/A29_Histogramme.csv","w");
    gsl_histogram_fprintf (file_histogram, hist_with_gsl, "%g", "%g");
    fclose(file_histogram);
    free(rng);
    return 0;
}