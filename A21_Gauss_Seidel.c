#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cvc_numerics.h"


// Physikalische Parameter des Systems
double BOX_X = 1;
double BOX_Y = 1;
double V_RAND = 0;

double R_ELECTRODE = 0.1;

double x_LEFT_ELECTRODE = 0.25;
double x_RIGHT_ELECTRODE = 0.75;
double y_LEFT_ELECTRODE = 0.5;
double y_RIGHT_ELECTRODE = 0.5;

double V_LEFT_ELECTRODE = -1;
double V_RIGHT_ELECTRODE = 1;


// Typ eines Gitterpunktes
enum Node_Type {
    INSIDE,
    DIRECHLET,
    NEUMANN,
    DISABLED
};


// Struktur eines Gitterpunktes: Typ und Wert
struct node {
    enum Node_Type type;
    double value;
};


// Initialisierung eines quadratischen Arrays F[N][N]
int InitDomain(int N, struct node F[][N], double N_val, double E_val, double S_val, double W_val) {
    for (int n_x = 0; n_x < N; n_x++) {
        for (int n_y = 0; n_y < N; n_y++) {
            if (n_x == 0 || n_x == N-1 || n_y == 0 || n_y == N-1) {
                F[n_x][n_y].type = DIRECHLET;
                F[n_x][n_y].value = V_RAND;
            } else {
                F[n_x][n_y].type = INSIDE;
                F[n_x][n_y].value = 0;
            }
        }
    }
    return 0;
}


// Einsetzten eines Potentialkreises mit Radius < R um c_x, c_y
int InsertCircle(int N, double c_x, double c_y, double R, struct node F[][N], double val) {
    double delta_x = BOX_X / N;
    double delta_y = BOX_Y / N;
    for (int n_x = 0; n_x < N; n_x++) {
        for (int n_y = 0; n_y < N; n_y++) {
            if (cvc_norm_2D(delta_x * n_x, delta_y * n_y) < R) {
                F[n_x][n_y].type = DIRECHLET;
                F[n_x][n_y].value = val;
            }
        }
    }
    return 0;
}


// Gauß-Seidel-Verfahren für die Poisson-Gleichung
void PoissonGaussSeidel(int N, struct node F[][N], int iter_max, double tolerance) {
    int count = 0;                                                                  // Anzahl der Schleifendurchläufe
    double R_squared = 0;                                                           // quadratische relative Anderung aller Gitterpunkte
    while (count < iter_max || sqrt(R_squared) > tolerance) {
        for (int n_x = 0; n_x < N; n_x++) {
            for (int n_y = 0; n_y < N; n_y++) {   
                switch (F[n_x][n_y].type) {
                    case INSIDE:
                        double F_old = F[n_x][n_y].value;
                        double F_new = 1 / 4 * (F[n_x][n_y-1].value + F[n_x][n_y+1].value + F[n_x-1][n_y].value + F[n_x+1][n_y].value);
                        F[n_x][n_y].value = F_new;
                        R_squared += cvc_npow((F_new-F_old) / F_old, 2);            // neue relative Änderung für innere Werte
                        break;
                    case DIRECHLET:
                        break;
                    default:
                        break;
                }
            }
        }
        count++;
    }

    // // Einssetzen des x-Vektors
    // double *x_vector = (double*) calloc(sizeof(double) * N);
    // for (int x_i = 0; x_i < N; x_i++) {
    //     x_vector[x_i] = 1;
    // }

    // for (int x_i = 0; x_i < N; x_i++) {
    //     for (L_i = 0; L_i < x_i; L_i++) {
    //         x_vector[x_i] += F[x_i][L_i] * x_vector[L_i];
    //     }
    //     for (R_i = x_i + 1; R_i < N; R_i++) {
    //         x_vector[x_i] += F[x_i][R_i] * x_vector[R_i];
    //     }
    //     x_vector[x_i] = 
    // }
    // return;
}

int main(void) {
    // Parameter des Gauß-Seidel-Verfahrens
    double tolerance = 1e-3;
    int iter_max = 1000;

    // Files
    FILE* field_file = fopen("data/A21_Gauss_Seidel_N32_data.csv", "w");

    int N = 32;
    struct node (*F)[N] = malloc(sizeof(struct node[N][N]));                                            // Erstellen des Feldes
    InitDomain(N, F, 0, 0, 0, 0);                                                                  // Initialisierung des Randes
    InsertCircle(N, x_LEFT_ELECTRODE, y_LEFT_ELECTRODE, R_ELECTRODE, F, V_LEFT_ELECTRODE);         // Einfügen linke Elektrode
    InsertCircle(N, x_RIGHT_ELECTRODE, y_RIGHT_ELECTRODE, R_ELECTRODE, F, V_RIGHT_ELECTRODE);      // Einfügen rechte Elektrode
    PoissonGaussSeidel(N, F, iter_max, tolerance);                                                 // Anwendung des Gauß-Seidel-Verfahrens

    for (int n_x = 0; n_x < N; n_x++) {
        for (int n_y = 0; n_y < N; n_y++) {   
            fprintf(field_file, "%g\n", F[n_x][n_y].value);
        }
    }
    fclose(field_file);
    free(F);


    // F[1][0].value = 5;
    // printf("%g\n", F[1][0].value);
    // for (int N = 32; N <= 256; N *= 2) {                    // Itieren über Gittergrößen N = 32, 64, 128 und 256

    // }

    return 0;
}