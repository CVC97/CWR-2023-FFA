#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cvc_numerics.h"


// Physikalische Parameter des Systems
const double BOX_X = 1;
const double BOX_Y = 1;
const double V_RAND = 0;


// Elektroden
const double R_ELECTRODE = 0.1;
const double x_LEFT_ELECTRODE = 0.25;          
const double x_RIGHT_ELECTRODE = 0.75;
const double y_LEFT_ELECTRODE = 0.5;
const double y_RIGHT_ELECTRODE = 0.5;
const double V_LEFT_ELECTRODE = -1;
const double V_RIGHT_ELECTRODE = 1;


// Hilfsparameter
const double SMOOTHING_FACTOR = 1e-10;


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
            if (cvc_norm_2D(delta_x * n_x - c_x, delta_y * n_y - c_y) < R) {
                F[n_x][n_y].type = DIRECHLET;
                F[n_x][n_y].value = val;
            }
        }
    }
    return 0;
}


// Gauß-Seidel-Verfahren für die Poisson-Gleichung
int PoissonGaussSeidel(int N, struct node F[][N], int iter_max, double tolerance) {
    int count = 0;                                                                                      // Anzahl der Schleifendurchläufe
    double R_squared;                                                                                   // quadratische relative Anderung aller Gitterpunkte
    do {
        R_squared = 0;
        for (int n_x = 0; n_x < N; n_x++) {
            for (int n_y = 0; n_y < N; n_y++) {   
                switch (F[n_x][n_y].type) {
                    case INSIDE:
                        double F_old = F[n_x][n_y].value;
                        double F_new = 1.0 / 4 * (F[n_x][n_y-1].value + F[n_x][n_y+1].value + F[n_x-1][n_y].value + F[n_x+1][n_y].value);
                        F[n_x][n_y].value = F_new;
                        R_squared += cvc_npow((F_new-F_old) / (F_old + SMOOTHING_FACTOR), 2);           // neue relative Änderung für innere Werte
                        break;

                    case DIRECHLET:
                        break;

                    default:
                        break;
                }
            }
        }
        count++;
    } while (count < iter_max && sqrt(R_squared) > tolerance);
    return count;
}


// SOR-Verfahren für die Poisson-Gleichung
int PoissonSOR(int N, struct node F[][N], int iter_max, double tolerance, double alpha) {
    int count = 0;                                                                                      // Anzahl der Schleifendurchläufe
    double R_squared;                                                                                   // quadratische relative Anderung aller Gitterpunkte
    do {
        R_squared = 0; 
        for (int n_x = 0; n_x < N; n_x++) {
            for (int n_y = 0; n_y < N; n_y++) {   
                switch (F[n_x][n_y].type) {
                    case INSIDE:
                        double F_old = F[n_x][n_y].value;
                        double F_new = 1.0 / 4 * (F[n_x][n_y-1].value + F[n_x][n_y+1].value + F[n_x-1][n_y].value + F[n_x+1][n_y].value);
                        F[n_x][n_y].value = F_old + alpha * (F_new - F_old);
                        R_squared += cvc_npow((F_new-F_old) / (F_old + SMOOTHING_FACTOR), 2);                                // neue relative Änderung für innere Werte
                        break;

                    case DIRECHLET:
                        break;

                    default:
                        break;
                }
            }
        }
        count++;
    } while (count < iter_max && sqrt(R_squared) > tolerance);
    return count;
}


int main(void) {
    // Parameter des SOR-Verfahrens
    double tolerance = 1e-4;
    int iter_max = 1e6;
    int N = 16;

    FILE* SOR_convergence_file = fopen("data/A24_SOR_N16_convergence_data.csv", "w");
    fprintf(SOR_convergence_file, "alpha, iterations till convergence SOR, iterations till convergence Gauss-Seidel\n");

    struct node (*F_SOR)[N] = malloc(sizeof(struct node[N][N]));                                                    // Erstellen des Feldes


    // Feld für SOR und Anzahl Iterationen bis Konvergenz
    for (double alpha = 1; alpha < 2 - 0.01; alpha += 0.025) {
        InitDomain(N, F_SOR, 0, 0, 0, 0);                                                                           // Initialisierung des Randes
        InsertCircle(N, x_LEFT_ELECTRODE, y_LEFT_ELECTRODE, R_ELECTRODE, F_SOR, V_LEFT_ELECTRODE);                  // Einfügen linke Elektrode
        InsertCircle(N, x_RIGHT_ELECTRODE, y_RIGHT_ELECTRODE, R_ELECTRODE, F_SOR, V_RIGHT_ELECTRODE);               // Einfügen rechte Elektrode
        int n_iterations_SOR = PoissonSOR(N, F_SOR, iter_max, tolerance, alpha);                                    // Anwendung des SOR-Verfahrens
        printf("N%d\talpha: %g\n", N, alpha);
        fprintf(SOR_convergence_file, "%g, %d\n", alpha, n_iterations_SOR);                                         // Konvergenz für alpha                  
    }

    // fclose(SOR_convergence_file);
    free(F_SOR);
    return 0;
}