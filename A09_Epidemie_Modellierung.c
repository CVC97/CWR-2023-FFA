#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cvc_numerics.h"


// Pandemie-Parameter
const double alpha = 1 / 6.1;
const double Epi_Mod_gamma = 1 / 3.7;

const double N = 83e6;
const double E0 = 30000;
const double I0 = 9000;


// SEIR-DGL
int SEIR_ODE(double t, const double y[], double f[], void *params) {
    // Parameter alpha, beta, gamma
    double SEIR_alpha = ((double*) params)[0];
    double SEIR_beta = ((double*) params)[1]; 
    double SEIR_gamma = ((double*) params)[2]; 

    // Berechnung der Gesamtpopulation aus SEIR
    int SEIR_N = 0;                                                     
    for (int i = 0; i < 4; i++) {
        SEIR_N += y[i];
    }

    // Berechnung der Ableitungen von SEIR
    f[0] = -SEIR_beta * y[2] * y[0] / SEIR_N;                           // dS/dt (Änderung Susceptible)
    f[1] = SEIR_beta * y[2] * y[0] / SEIR_N - SEIR_alpha * y[1];        // dE/dt (Änderung Exposed)
    f[2] = SEIR_alpha * y[1] - SEIR_gamma * y[2];                       // dI/dt (Änderung Infected)             
    f[3] = - (f[0] + f[1] + f[2]);                                      // Änderung Removed
    return 0;
}


int main(void) {
    double S = N - (E0 + I0), E = E0, I = I0, R = 0, t = 0, I_max_1_25 = I0, I_max_1_5 = I0, I_max_2 = I0;

    // R0 = 1.25
    double params_SEIR_1_25[3] = {alpha, 1.25*Epi_Mod_gamma, Epi_Mod_gamma};
    double y_SEIR_1_25[4] = {S, E, I, R};

    FILE* file_epidemic_modelling_1_25 = fopen("data/A09_Epidemie_Modellierung_R0_1.25.csv", "w");
    fprintf(file_epidemic_modelling_1_25, "Zeit t in Tagen, Susceptible S, Exposed E, Infectious I, Removed R\n");
    fprintf(file_epidemic_modelling_1_25, "%f, %f, %f, %f, %f\n", t, y_SEIR_1_25[0], y_SEIR_1_25[1], y_SEIR_1_25[2], y_SEIR_1_25[3]);

    // R0 = 1.5
    double params_SEIR_1_5[3] = {alpha, 1.5*Epi_Mod_gamma, Epi_Mod_gamma};
    double y_SEIR_1_5[4] = {S, E, I, R};

    FILE* file_epidemic_modelling_1_5 = fopen("data/A09_Epidemie_Modellierung_R0_1.5.csv", "w");
    fprintf(file_epidemic_modelling_1_5, "Zeit t in Tagen, Susceptible S, Exposed E, Infectious I, Removed R\n");
    fprintf(file_epidemic_modelling_1_5, "%f, %f, %f, %f, %f\n", t, y_SEIR_1_5[0], y_SEIR_1_5[1], y_SEIR_1_5[2], y_SEIR_1_5[3]);

    // R0 = 2
    double params_SEIR_2[3] = {alpha, 2*Epi_Mod_gamma, Epi_Mod_gamma};
    double y_SEIR_2[4] = {S, E, I, R};

    FILE* file_epidemic_modelling_2 = fopen("data/A09_Epidemie_Modellierung_R0_2.csv", "w");
    fprintf(file_epidemic_modelling_2, "Zeit t in Tagen, Susceptible S, Exposed E, Infectious I, Removed R\n");
    fprintf(file_epidemic_modelling_2, "%f, %f, %f, %f, %f\n", t, y_SEIR_2[0], y_SEIR_2[1], y_SEIR_2[2], y_SEIR_2[3]);


    // Integrations-Parameter
    double delta_t = 0.01;
    while (t < 365) {
        t += delta_t;
        cvc_euler_step(t, delta_t, y_SEIR_1_25, SEIR_ODE, 4, params_SEIR_1_25);
        cvc_euler_step(t, delta_t, y_SEIR_1_5, SEIR_ODE, 4, params_SEIR_1_5);
        cvc_euler_step(t, delta_t, y_SEIR_2, SEIR_ODE, 4, params_SEIR_2);
        fprintf(file_epidemic_modelling_1_25, "%f, %f, %f, %f, %f\n", t, y_SEIR_1_25[0], y_SEIR_1_25[1], y_SEIR_1_25[2], y_SEIR_1_25[3]);
        fprintf(file_epidemic_modelling_1_5, "%f, %f, %f, %f, %f\n", t, y_SEIR_1_5[0], y_SEIR_1_5[1], y_SEIR_1_5[2], y_SEIR_1_5[3]);
        fprintf(file_epidemic_modelling_2, "%f, %f, %f, %f, %f\n", t, y_SEIR_2[0], y_SEIR_2[1], y_SEIR_2[2], y_SEIR_2[3]);
        if (y_SEIR_1_25[2] > I_max_1_25) {
            I_max_1_25 = y_SEIR_1_25[2];
        }
        if (y_SEIR_1_5[2] > I_max_1_5) {
            I_max_1_5 = y_SEIR_1_5[2];
        }
        if (y_SEIR_2[2] > I_max_2) {
            I_max_2 = y_SEIR_2[2];
        }
    }
    printf("Maximale Anzahl Infizierter (R0 = 1.25): %f\n", I_max_1_25);
    printf("Maximale Anzahl Infizierter (R0 = 1.5): %f\n", I_max_1_5);
    printf("Maximale Anzahl Infizierter (R0 = 2): %f\n", I_max_2);
    fclose(file_epidemic_modelling_1_25);
    fclose(file_epidemic_modelling_1_5);
    fclose(file_epidemic_modelling_2);
    return 0;
}