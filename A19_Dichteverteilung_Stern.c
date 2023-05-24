#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cvc_numerics.h"


// Physikalische Parameter des Systems
const double GAMMA = 3.0 / 4;
const double KAPPA = 3.85e9;
const double R_SUN = 7e8;
const double G = 6.67384e-11;


// Hilfsparameter
const double DIMENSION = 2;


// Berechnung der Masse- und Dichteänderung für gegebene Radii
int ODE_star_density(double r, const double y[], double f[], void *params) {
    double rho_0 = ((double*) params)[0];
    double epsilon = ((double*) params)[1];

    // Berechnung der Masse des Sterns von 0 bis r
    double inner_mass = 4.0 / 3 * cvc_PI * rho_0 * cvc_npow(epsilon, 3);                    // Masse von 0 bis epsilon (Dichte konstant)
    double outer_mass = 0;                                                                  // Masse von epsilon bis r mit RK4
    double dr = 1e3;                                                                        // dr für die Masseintegration
    while (epsilon + dr < r) {

        epsilon += dr;
    }
    if (r < epsilon) {
        f[0] = 0:
    } else {
        f[0] = 3.0 / 2 * 1 / (KAPPA * GAMMA) * pow(rho_0, 2-GAMMA) * G * m / cvc_npow(r, 2);
    }

    return 0;
}


int main(void) {
    double R = 1.2 * R_SUN;
    double rho_0 = 7e4;
    double epsilon = 1e-4;

    double params[2] = {rho_0, epsilon};


    return 0;
}