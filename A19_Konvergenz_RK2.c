#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cvc_numerics.h"


// physikalische Parameter und Konstanten
const int N = 1;
const double pi = 3.14159265359;
const double k = 100;
const double mass = 1;
const double v_0 = 1;


int F_pendulum_ode(double t, const double y[], double f[], void *params) {
    f[0] = y[1];
    f[1] = -k*y[0] / mass;
    return 0;
}


int main(void) {
    int dimension = 2 * N;

    // Feste analytische Parameter
    double T = 2 * pi * sqrt(mass / k);
    double x_T4 = v_0 * sqrt(mass / k);
    double T4 = T/4;

    // Zustandsarray
    double y_rk2[2] = {0, v_0};
    double y_rk4[2] = {0, v_0};
    double y_euler[2] = {0, v_0};
    double y_verlet[2] = {0, v_0};

    printf("x_T4: %g, T4: %g\n", x_T4, T4);

    // Erstellen der Residuendatei
    FILE* res_file = fopen("data/A19_Konvergenz_RK2.csv", "w");
    fprintf(res_file, "Zeitschritt delta_t, Residuum RK2, , Residuum RK4, Residuum Euler, Residuum Verlet\n");

    // Integrationsparameter
    double t = 0, d_t, exp = 0, residue_RK2, residue_RK4, residue_Euler, residue_verlet;
    for (int i = 0; i < 100; i++) {
        t = 0;                                                          // Starten der jeweiligen Integration bei t = 0
        exp = -1 - i * (7.0/99);                                        // Anpassen des Exponentens des Schrittweite
        d_t = pow(10, exp);                                             // neue Schrittweite aus Exponent

        // Numerische Integration bis vorletztem Integrationsschritt
        while (t + d_t < T4) {
            rk2_step(t, d_t, y_rk2, F_pendulum_ode, dimension, NULL);
            //printf("x (RK29: %g\n", y_rk2[0]);
            rk4_step(t, d_t, y_rk4, F_pendulum_ode, dimension, NULL);
            euler_step(t, d_t, y_euler, F_pendulum_ode, dimension, NULL);
            verlet_step(t, d_t, y_verlet, F_pendulum_ode, dimension, NULL);
            t += d_t;
        }

        // letzter Integrationsschritt mit verbleibender Zeit
        rk2_step(t, T4-t, y_rk2, F_pendulum_ode, dimension, NULL);
        rk2_step(t, T4-t, y_rk4, F_pendulum_ode, dimension, NULL);
        euler_step(t, T4-t, y_euler, F_pendulum_ode, dimension, NULL);
        verlet_step(t, T4-t, y_verlet, F_pendulum_ode, dimension, NULL);

        // Berechnung und Ausgabe des "Residuums"         
        residue_RK2 = fabs(y_rk2[0] - x_T4);
        printf("Resi_RK2: %g\n", fabs(y_rk2[0] - x_T4));
        residue_RK4 = fabs(y_rk4[0] - x_T4);
        residue_Euler = fabs(y_euler[0] - x_T4);
        residue_verlet = fabs(y_verlet[0] - x_T4);
        fprintf(res_file, "%g, %g, %g, %g, %g\n", d_t, residue_RK2, residue_RK4, residue_Euler, residue_verlet);
    }
    
    fclose(res_file);
    return 0;
}