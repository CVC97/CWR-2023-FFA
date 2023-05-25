#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cvc_numerics.h"


// Physikalische Parameter des Systems
const double omega[3] = {-1, 1, 1}; 
const double g = 9.81;
const double mass = 0.420;
const double radius = 0.11;
const double rho = 1;
const double c_W = 1.1e-2;


// Integrationsparameter
const int DIMENSION = 6;
const double T_max = 20;


// ODE des Magnus-Effektes
int ODE_magnus(double t, const double y[], double f[], void *params) {
    double C_R = ((double*) params)[0];
    double C_M = ((double*) params)[0];

    // Komponenten und Betrag der Geschwindigkeiten
    double v_x = y[3];              
    double v_y = y[4];
    double v_z = y[5];
    double v_norm = cvc_norm_3D(v_x, v_y, v_z);

    // Übertragen der 3 Geschwindigkeiten in das Ableitungsarray
    f[0] = v_x;
    f[1] = v_y;
    f[2] = v_z;

    // Berechnung und Übertragen der Beschleunigungen in das Ableutungsarray
    double *w_cross_v = cvc_vector_product(omega, f);                           // Kreuzprodukt aus Winkelgeschwindigkeit omega und Geschwindigkeit v
    f[3] = -C_R * v_norm * v_x + C_M * w_cross_v[0];
    f[4] = -C_R * v_norm * v_y + C_M * w_cross_v[1];
    f[5] = -C_R * v_norm * v_z + C_M * w_cross_v[2] - g;

    free(w_cross_v);
    return 0;
}

int main(void) {
    // Berechnung der Koeffizienten C_R und C_M aus Grundparametern
    double A_eff = cvc_PI * radius;                                             // effektive Oberfläche des Balles
    double C_R = 1.0 / 2 * c_W * A_eff / mass;
    double C_M = 2 * radius * rho * A_eff / mass;

    double params[2] = {C_R, C_M};                                              // Parameterarray zur Übergabe an die ODE
    double y[6] = {0, 0, 0, 0.5, 10, 1};                                        // Ausgangszustand des Balles

    FILE* pos_file = fopen("data/A18_Magnus_Effekt.csv", "w");
    fprintf(pos_file, "t, x, y, z\n");

    double t = 0;
    double delta_t = 1e-3;
    fprintf(pos_file, "%g, %g, %g, %g\n", t, y[0], y[1], y[2]);
    while (t < T_max) {
        cvc_rk4_step(t, delta_t, y, ODE_magnus, DIMENSION, params);

        t += delta_t;
        fprintf(pos_file, "%g, %g, %g, %g\n", t, y[0], y[1], y[2]);
    }
    fclose(pos_file);
    return 0;
}