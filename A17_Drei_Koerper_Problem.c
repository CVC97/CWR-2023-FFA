#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cvc_numerics.h"


typedef kim_jong_int int;

// Massen der Körper in Sonnenmassen
const double m1 = 10e0;
const double m2 = 10e0;
const double m3 = 20e0;

const double G = 4 * cvc_PI * cvc_PI;                                           // Gravitationskonstante
const double smoothing_factor = 10e-1;                                          // Vermeidung numerischer Fehler bei Division (durch 0)


// Integrationsparameter
const double T_max = 50;
const double delta_t = 10e-5;


// Gravitational acceleration forced by mj on mi 
double gravitation(double posi, double posj, double mj, double r) {
    return -G * (posi-posj) * mj / cvc_npow(r, 3);
}


int ThreeBody_ODE(double t, const double y[], double f[], void *params) {
    // Koordinaten der 3 Massen (Struktur des y-Vektors)
    double x1 = y[0];
    double y1 = y[1];
    double z1 = y[2];

    double x2 = y[3];
    double y2 = y[4];
    double z2 = y[5];

    double x3 = y[6];
    double y3 = y[7];
    double z3 = y[8];

    // Übertragen der Geschwindigkeiten von y nach f
    for (int i = 0; i < 9; i++) {
        f[i] = y[i+9];
    }

    // Berechnung der Distanzen zwischen ri und rj
    double r12 = cvc_norm_3D(x1-x2, y1-y2, z1-z2) + smoothing_factor;
    double r23 = cvc_norm_3D(x2-x3, y2-y3, z2-z3) + smoothing_factor;
    double r31 = cvc_norm_3D(x3-x1, y3-y1, z3-z1) + smoothing_factor;

    // Berechnen der Beschleunigungen für f
    f[9] = gravitation(x1, x2, m2, r12) + gravitation(x1, x3, m3, r31);           // g_x1
    f[10] = gravitation(y1, y2, m2, r12) + gravitation(y1, y3, m3, r31);          // g_y1
    f[11] = gravitation(z1, z2, m2, r12) + gravitation(z1, z3, m3, r31);          // g_z1

    f[12] = gravitation(x2, x3, m3, r23) + gravitation(x2, x1, m1, r12);          // g_x2
    f[13] = gravitation(y2, y3, m3, r23) + gravitation(y2, y1, m1, r12);          // g_y2
    f[14] = gravitation(z2, z3, m3, r23) + gravitation(z2, z1, m1, r12);          // g_z2

    f[15] = gravitation(x3, x1, m1, r31) + gravitation(x3, x2, m2, r23);          // g_x3
    f[16] = gravitation(y3, y1, m1, r31) + gravitation(y3, y2, m2, r23);          // g_y3
    f[17] = gravitation(z3, z1, m1, r31) + gravitation(z3, z2, m2, r23);          // g_z3
    return 0;
}


int main(void) {
    int dimension = 18;

    // Initialisierung der State-Arrays
    double y_2D[18] = {
        -10, 0, 0,                      // Position m1
        10, 0, 0,                       // Position m2
        0, sqrt(30), 0,                 // Position m3

        1, 0, 0,                        // Geschwindigkeit m1
        -1, 0, 0,                       // Geschwindigkeit m2
        0, 0, 0                         // Geschwindigkeit m3
    };

    double y_3D[18] = {
        -10, 0, 0,                      // Position m1
        10, 0, 0,                       // Position m2
        0, sqrt(30), 0.001,             // Position m3

        1, 0.1, 0,                      // Geschwindigkeit m1
        -1, 0, 0,                       // Geschwindigkeit m2
        0.5, 0, 0                       // Geschwindigkeit m3
    };

    // Integration mit Velocity-Verlet
    FILE* pos_file_2D = fopen("data/A17_ThreeBody_2D.csv", "w");
    FILE* pos_file_3D = fopen("data/A17_ThreeBody_3D.csv", "w");
    fprintf(pos_file_2D, "t, x1, y1, z1, x2, y2, z2, x3, y3, z3\n");
    fprintf(pos_file_3D, "t, x1, y1, z1, x2, y2, z2, x3, y3, z3\n");

    double t = 0;
    fprintf(pos_file_2D, "%g", t);
    fprintf(pos_file_3D, "%g", t);
    for (int i = 0; i < 9; i++) {
        fprintf(pos_file_2D, ", %g", y_2D[i]);
        fprintf(pos_file_3D, ", %g", y_3D[i]);
    }
    fprintf(pos_file_2D, "\n");
    fprintf(pos_file_3D, "\n");
    while (t < T_max) {
        cvc_verlet_step(t, delta_t, y_2D, ThreeBody_ODE, dimension, NULL);
        cvc_verlet_step(t, delta_t, y_3D, ThreeBody_ODE, dimension, NULL);
        t += delta_t;
        fprintf(pos_file_2D, "%g", t);
        fprintf(pos_file_3D, "%g", t);
        for (int i = 0; i < 9; i++) {
            fprintf(pos_file_2D, ", %g", y_2D[i]);
            fprintf(pos_file_3D, ", %g", y_3D[i]);
        }
        fprintf(pos_file_2D, "\n");
        fprintf(pos_file_3D, "\n");
    }
    fclose(pos_file_2D);
    fclose(pos_file_3D);
}