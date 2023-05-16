#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cvc_numerics.h"


// Physiklische Konstanten
const double epsilon_0 = 8-85419e-12;
const double e_positron = 1.602e-19;
const double m_positron = 9.109e-31;


// Coublombkraft für gegebenen Abstand r (Radialsymmetrie)
int ODE_Coulomb(double t, const double y[], double f[], void *params) {
    double F_Coulomb = cvc_npow(e_positron, 2) / (4*cvc_PI*epsilon_0) / cvc_npow(y[0], 2);
    f[0] = y[1];
    f[1] = F_Coulomb / m_positron;
    return 0;
}


// Lösung des Anfangswertproblems mit Verlet, Rückgabe Differenz zur analytischen r_tf
double G(double v0) {
    // Zeit- und Ortsparameter
    double t0 = 0;                                          // Startzeit t0
    double tf = 2;                                          // Endzeit tf
    double r0 = 10;                                         // Startposition r_t0
    double rf = 10;                                         // Endposition r_tf

    // Numerische Integration mit Verlet
    double y[2] = {r0, v0};                                 // Initialisierung des Startarrays
    double t = 0, delta_t = 10e-7;                          // Zeitparameter (geeignetes delta_t nach Konvergenz des Verlet)
    while (t + delta_t < tf) {
        cvc_verlet_step(t, delta_t, y, ODE_Coulomb, 2, NULL);
        t += delta_t;
    }
    cvc_verlet_step(t, tf - t, y, ODE_Coulomb, 2, NULL);    // Ausgleichszeitschritt
    return (y[0] - rf);
}


int main(void) {
    double v0 = -5;                                         // Startgeschwindigkeit v_t0

    // Parameter der Newton-Raphson
    double delta = 10e-3;                                   // verwendetes Delta bei der Ableitung der Funktion
    double rel_tol = 10e-2;                                 // relative Fehlertoleranz
    int max_iter = 10e3;                                    // maximale Anzahl an Iterationen im NR

    // Aufruf der Newton-Raphson
    double root_G = cvc_find_root_newton_raphson(G, v0, delta, rel_tol, max_iter);
    printf("Nullstelle von G(v0): v0 = %g\n", root_G);

    // Variation des v0 in [-12, -2]
    FILE* diff_file = fopen("data/old_A23_Shooting_Methode_Gv.csv", "w");
    fprintf(diff_file, "v0, delta\n");
    int n_v_steps = 100;
    for (int i = 0; i < n_v_steps; i++) {
        v0 = -2 - 10.0*i/(n_v_steps - 1);
        printf("G(v0) = %g\t v0 = %g\n", G(v0), v0);
        fprintf(diff_file, "%g, %g\n", v0, G(v0));
    }
    fclose(diff_file);

    // Berechnung der Trajektorien mit Verlet
    return 0;
}