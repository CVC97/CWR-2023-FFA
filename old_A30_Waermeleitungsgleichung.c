#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cvc_numerics.h"


const double L = 1;                                     // Länge der Domäne
const double delta_x = 0.01;                            // Diskretisierung der Domäne
const double D = 0.1;                                   // Wärmediffusivität
const int N = L / delta_x;                              // Anzahl der Diskretisierungsschritte


int heat_FTCS(double t, double y[], double dt) {
    double *y_temp = (double*) malloc(N * sizeof(double));
    for(size_t i = 0; i < N; ++i) {
        double T = y[i % N];
        double TW = y[(i - 1) % N];
        double TO = y[(i + 1) % N];

        y_temp[i % N] = 1;
    }

    return 0;
}

int main(void) {

    return 0;
}