#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cvc_numerics.h"


// physikalische Parameter und Konstanten
const int N = 1;
const double pi = 3.14159265359;
const double k = 100;
const double mass = 1;
const double g = 9.81;
const double v_0 = 1;


int F_pendulum_ode(double t, const double y[], double f[], void *params) {
    // Reibungsparameter
    double v = y[1];
    double mu_s = ((double*) params)[0];
    double mu_k = ((double*) params)[1];
    double b = ((double*) params)[2];  

    // Berechnung der Kräfte
    double F = -k*y[0];
    double F_Normal = mass * g;
    double F_Haft = -mu_s * F_Normal;
    double F_Gleit = -mu_k * F_Normal * v / fabs(v);
    double F_viskos = -b * v;

    // Berechnung des Ableitungsvektors
    f[0] = v;
    if (fabs(v) < 0.01) {                    
        if (fabs(F_Haft) > fabs(F + F_Gleit + F_viskos)) {                // Keine Beschleunigung wenn Haftreibung größer als Summer restlicher Reibungen
            f[1] = 0;
            return 1;
        }
    }
    f[1] = (F + F_Gleit + F_viskos) / mass;
    return 0;
}


int main(void) {
    int dimension = 2*N;
    double omega_0 = sqrt(k/mass);
    double t = 0, delta_t = 1e-4, T_max = 2.5;

    double params_no_friction[3] = {0, 0, 0};
    double y_no_friction[2] = {0, v_0};

    // Haftreibung
    double params_hafti_a[3] = {1, 0, 0};
    double params_hafti_b[3] = {2, 0, 0};

    double y_hafti_a[2] = {0, v_0};
    double y_hafti_b[2] = {0, v_0};


    // Gleitreibung
    double params_gleit_a[3] = {0, 0.1, 0};
    double params_gleit_b[3] = {0, 0.2, 0};

    double y_gleit_a[2] = {0, v_0};
    double y_gleit_b[2] = {0, v_0};


    // Viskose Reibung
    double params_viskos_a[3] = {0, 0, 2*mass*omega_0 / 4};
    double params_viskos_b[3] = {0, 0, 2*mass*omega_0};
    double params_viskos_c[3] = {0, 0, 2*mass*omega_0 * 4};

    double y_viskos_a[2] = {0, v_0};
    double y_viskos_b[2] = {0, v_0};
    double y_viskos_c[2] = {0, v_0};

    FILE* pos_file = fopen("data/old_A20_Reibung.csv", "w");
    fprintf(pos_file, "Zeit t, keine Reibung, Gleit a, Gleit b, Viskos a, Viskos b, Viskos c\n");
    fprintf(pos_file, "%g, %g, %g, %g, %g, %g, %g\n", t, y_no_friction[0], y_gleit_a[0], y_gleit_b[0], y_viskos_a[0], y_viskos_b[0], y_viskos_c[0]);


    // Integration mit RK2
    while (t < T_max) {
        rk2_step(t, delta_t, y_no_friction, F_pendulum_ode, dimension, params_no_friction);

        rk2_step(t, delta_t, y_gleit_a, F_pendulum_ode, dimension, params_gleit_a);
        rk2_step(t, delta_t, y_gleit_b, F_pendulum_ode, dimension, params_gleit_b);

        rk2_step(t, delta_t, y_viskos_a, F_pendulum_ode, dimension, params_viskos_a);
        rk2_step(t, delta_t, y_viskos_b, F_pendulum_ode, dimension, params_viskos_b);
        rk2_step(t, delta_t, y_viskos_c, F_pendulum_ode, dimension, params_viskos_c);

        fprintf(pos_file, "%g, %g, %g, %g, %g, %g, %g\n", t, y_no_friction[0], y_gleit_a[0], y_gleit_b[0], y_viskos_a[0], y_viskos_b[0], y_viskos_c[0]);
        t += delta_t;
    }

    fclose(pos_file);
    return 0;
}