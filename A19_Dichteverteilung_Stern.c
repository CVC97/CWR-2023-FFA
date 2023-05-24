#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cvc_numerics.h"


// Physikalische Parameter des Systems
const double GAMMA = 4.0 / 3;
const double KAPPA = 3.85e9;
const double R_SUN = 7e8;
const double G = 6.67384e-11;

// Chemische Paramater
const double UNIVERSAL_GAS_CONSTANT = 8.314;
const double VOLUME_MOL = 22.4e-3;


// Integrationsparameter
const double EPSILON = 1e-4;
const int DIMENSION = 2;


// Berechnung der Masse- und Dichteänderung für gegebenen Radius, sowie Masse-, und Dichtezustand
int PDE_star_density(double r, const double y[], double f[], void *params) {
    double rho = y[0];                                                                              // Dichte
    double mass = y[1];                                                                             // Masse

    f[0] = -3.0 / 2 / (KAPPA * GAMMA) * pow(rho, 2-GAMMA) * G * mass / cvc_npow(r, 2);              // Dichteänderung
    f[1] = 4 * cvc_PI * rho * cvc_npow(r, 2);                                                       // Masseänderung
    return 0;
}


// Differenz numerischer und numerischer Radius (bei Erreichen der Grenzdichte) für gegebene Zentrumsdichte rho_0
double Delta_R(double rho_0) {
    double R_analytical = 1.2 * R_SUN;                                                              // analytischer Radius des Sterns

    // Startbedingung der numerischen Integration
    double m_0 = 4.0 / 3 * cvc_PI * rho_0 * cvc_npow(EPSILON, 3);                                   // Startmasse m_0 aus rho_0 und EPSILON
    double y[2] = {rho_0, m_0};                                                                     // Array mit Ausgangszustand (rho_0, m_0) im (numerischen Zentrum)

    double r = EPSILON;
    double delta_r = 1000;
    while (y[0] > 0.1) {
        cvc_rk4_step(r, delta_r, y, PDE_star_density, DIMENSION, NULL);
        r += delta_r;
    }
    printf("r (numerisch): %g\n", r);
    printf("r (analytisch): %g\n", R_analytical);
    return fabs(R_analytical - r);
}


int main(void) {
    // Suchen der Nullstelle mit Newton_Raphson
    double rho_0 = 7e4;                                                                             // Startdichte zur Nullstellensuche
    double delta_rho = 1e3;                                                                         // delta_rho der numerischen Differenzierung in der Nullsellensuche
    double rel_tol = 1e-3;                                                                          // Relative Toleranz zum Abbruch
    double max_iter = 1e4;                                                                          // Maximale Anzahl der Schleifendurchläufe der Nullstellensuche

    double rho_0_numerical = cvc_find_root_newton_raphson(Delta_R, rho_0, delta_rho, rel_tol, max_iter);
    printf("Numerisches rho_0: %g\n", rho_0_numerical);

    printf("Delta_R: %g\n", Delta_R(rho_0));


    // Visuelle Bestimmung der Nullstelle


    // Berechnung der Dichteverteilung mit oben bestimmten rho_0_numerical
    FILE* density_file = fopen("data/A19_Dichteverteilung.csv", "w");
    fprintf(density_file, "r, density, mass, pressure, temperature\n");

    double R_analytical = 1.2 * R_SUN;                                                              // analytischer Radius des Sterns

    // Startbedingung der numerischen Integration
    double rho_0_correct = 6.7e4;
    double m_0 = 4.0 / 3 * cvc_PI * rho_0 * cvc_npow(EPSILON, 3);                                   // Startmasse m_0 aus rho_0 und EPSILON
    double y[2] = {rho_0_correct, m_0};                                                             // Array mit Ausgangszustand (rho_0, m_0) im (numerischen Zentrum)

    double r = EPSILON;
    double delta_r = 1000;

    double pressure = KAPPA * pow(y[0], GAMMA);
    double temperature = pressure / UNIVERSAL_GAS_CONSTANT * VOLUME_MOL;
    fprintf(density_file, "%g, %g, %g, %g, %g\n", r, y[0], y[1], pressure, temperature);
    while (y[0] > 0.1) {
        cvc_rk4_step(r, delta_r, y, PDE_star_density, DIMENSION, NULL);
        r += delta_r;

        pressure = KAPPA * pow(y[0], GAMMA);
        temperature = pressure / UNIVERSAL_GAS_CONSTANT * VOLUME_MOL;     
        fprintf(density_file, "%g, %g, %g, %g, %g\n", r, y[0], y[1], pressure, temperature);
    }
    return 0;
}