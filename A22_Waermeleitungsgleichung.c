#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cvc_numerics.h"


const double L = 1;                                     // Länge der Domäne
const double delta_x = 0.01;                            // Diskretisierung der Domäne
const double D = 0.1;                                   // Wärmediffusivität
const double T_LINKS = 1;                               // Direchlet-Randbedingung links: T = 1
const double T_RECHTS = -1;                             // Direchlet-Randbedingung rechts: T = -1
const int N = L / delta_x;                              // Anzahl der Diskretisierungsschritte


int heat_FTCS(double t, double y[], double dt) {
    double *y_temp = (double*) malloc(N * sizeof(double));
    for(size_t i = 0; i < N; i++) {
        double T = y[i % N];
        double TW = y[(i - 1) % N];
        double TO = y[(i + 1) % N];

        // Ränder
        y_temp[i] = y[i] + dt * D * (TW - 2*T + TO) / cvc_npow(delta_x, 2);
    }

    // Übertragen der neuen Werte in das Ursprungsaray
    for(int i = 0; i < N; i++) {
        y[i] = y_temp[i];
    }

    free(y_temp);
    return 0;
}

int main(void) {

    return 0;
}