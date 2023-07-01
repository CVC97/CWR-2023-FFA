#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "cvc_numerics.h"
#include "cvc_rng.h"


int N = 200;                            // Größe des Ensembles


// Stochastische DGL des Aktienmarktes
int SDE_stock(double t, const double y[], double f[], double g[], void* params) {
    struct cvc_tuple_2 *T = (struct cvc_tuple_2*) params;
    double mu = T->x1;
    double sigma = T->x2;

    // Berechnung der f und g-Arrays aus y und den obigen Parametern
    for (int i = 0; i < N; i++) {
        f[i] = mu * y[i];
        for (int j = 0; j < N; j++) {
            if (i == j) {
                g[i * N + j] = sigma * y[j];    
            } else {
                g[i * N + j] = sigma * y[j] * 0.05;
            }
        }
    }
    return 0;
}


int main(void) {
    // Paramter des Systems
    double K0 = 10;                     // Startwert
    double mu = 0.1;                    // erwarete jährliche Rendite
    double sigma = 0.1;                 // Volatilität
    double delta_t = 10e-3;             // Zeitschritt in Jahren
    double t = 0;                       // Startzeit
    double t_month = 1.0/12;            // Monatsperiode
    double T_max = 20;                  // Endzeit

    // 2er-Tupel der Parameter mu und sigma
    struct cvc_tuple_2 parameters;
    parameters.x1 = mu;
    parameters.x2 = sigma;

    // Datenfile
    FILE* stock_file = fopen("data/A33_Aktien_Ensemble.csv", "w");
    FILE* stat_file = fopen("data/A33_Aktien_Ensemble_Stats.csv", "w");
    fprintf(stock_file, "%g", t);

    // Initialisieren und Ausgabe der Startwerte des Ensebles
    double *K = (double*) malloc(sizeof(double) * N);
    for (int i = 0; i < N; i++) {
        K[i] = K0;
        fprintf(stock_file, ", %g", K[i]);
    }

    int counter = 0;
    while (t < T_max) {
        t += delta_t;
        if (t - t_month * counter + delta_t > t_month) {
            cvc_eulerMaruyama_step(t, t_month * (counter+1) - t, K, SDE_stock, N, &parameters);
            double K_mean = cvc_mean(K, N);
            double K_sigma = cvc_sigma(K, N);
            fprintf(stock_file, "\n%g", t);
            for (int i = 0; i < N; i++) {
                fprintf(stock_file, ", %g", K[i]);
            }
            fprintf(stat_file, "%g ,%g, %g\n", t, K_mean, K_sigma);
            counter++;
        } else {
            cvc_eulerMaruyama_step(t, delta_t, K, SDE_stock, N, &parameters);
            fprintf(stock_file, "\n%g", t);
            for (int i = 0; i < N; i++) {
                fprintf(stock_file, ", %g", K[i]);
            }
        }
    }
    free(K);
    fclose(stock_file);
    fclose(stat_file);
    return 0;
}