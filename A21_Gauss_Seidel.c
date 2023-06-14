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
void PoissonGaussSeidel(int N, struct node F[][N], int iter_max, double tolerance) {
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
}


int main(void) {
    // Parameter des Gauß-Seidel-Verfahrens
    double tolerance = 1e-3;
    int iter_max = 1000;

    // Files
    FILE* field_32_file = fopen("data/A21_Gauss_Seidel_N32_data.csv", "w");
    FILE* field_64_file = fopen("data/A21_Gauss_Seidel_N64_data.csv", "w");
    FILE* field_128_file = fopen("data/A21_Gauss_Seidel_N128_data.csv", "w");
    FILE* field_256_file = fopen("data/A21_Gauss_Seidel_N256_data.csv", "w");

    // Göße der Felder
    int N32 = 32;
    int N64 = 64;
    int N128 = 128;
    int N256 = 256;

    // Erstellen der Felder
    struct node (*F32)[N32] = malloc(sizeof(struct node[N32][N32]));
    struct node (*F64)[N64] = malloc(sizeof(struct node[N64][N64]));
    struct node (*F128)[N128] = malloc(sizeof(struct node[N128][N128]));
    struct node (*F256)[N256] = malloc(sizeof(struct node[N256][N256]));

    // Initialisierung der Ränder
    InitDomain(N32, F32, 0, 0, 0, 0);
    InitDomain(N64, F64, 0, 0, 0, 0);
    InitDomain(N128, F128, 0, 0, 0, 0);
    InitDomain(N256, F256, 0, 0, 0, 0);

    // Hinzufügen der Elektroden
    InsertCircle(N32, x_LEFT_ELECTRODE, y_LEFT_ELECTRODE, R_ELECTRODE, F32, V_LEFT_ELECTRODE);                      // Einfügen linke Elektrode
    InsertCircle(N32, x_RIGHT_ELECTRODE, y_RIGHT_ELECTRODE, R_ELECTRODE, F32, V_RIGHT_ELECTRODE);                   // Einfügen rechte Elektrode

    InsertCircle(N64, x_LEFT_ELECTRODE, y_LEFT_ELECTRODE, R_ELECTRODE, F64, V_LEFT_ELECTRODE);                      // Einfügen linke Elektrode
    InsertCircle(N64, x_RIGHT_ELECTRODE, y_RIGHT_ELECTRODE, R_ELECTRODE, F64, V_RIGHT_ELECTRODE);                   // Einfügen rechte Elektrode

    InsertCircle(N128, x_LEFT_ELECTRODE, y_LEFT_ELECTRODE, R_ELECTRODE, F128, V_LEFT_ELECTRODE);                    // Einfügen linke Elektrode
    InsertCircle(N128, x_RIGHT_ELECTRODE, y_RIGHT_ELECTRODE, R_ELECTRODE, F128, V_RIGHT_ELECTRODE);                 // Einfügen rechte Elektrode

    InsertCircle(N256, x_LEFT_ELECTRODE, y_LEFT_ELECTRODE, R_ELECTRODE, F256, V_LEFT_ELECTRODE);                    // Einfügen linke Elektrode
    InsertCircle(N256, x_RIGHT_ELECTRODE, y_RIGHT_ELECTRODE, R_ELECTRODE, F256, V_RIGHT_ELECTRODE);                 // Einfügen rechte Elektrode

    // Anwendung des Gauß-Seidel-Verfahrens
    printf("F32 calculating...\n");
    PoissonGaussSeidel(N32, F32, iter_max, tolerance);   
    printf("F64 calculating...\n");
    PoissonGaussSeidel(N64, F64, iter_max, tolerance);
    printf("F128 calculating...\n");
    PoissonGaussSeidel(N128, F128, iter_max, tolerance);
    printf("F256 calculating...\n");
    PoissonGaussSeidel(N256, F256, iter_max, tolerance);                                                   

    // Beschreiben der Files
    for (int n_x = 0; n_x < N32; n_x++) {
        for (int n_y = 0; n_y < N32; n_y++) {   
            fprintf(field_32_file, "%g\n", F32[n_x][n_y].value);
        }
    }
    for (int n_x = 0; n_x < N64; n_x++) {
        for (int n_y = 0; n_y < N64; n_y++) {   
            fprintf(field_64_file, "%g\n", F64[n_x][n_y].value);
        }
    }
    for (int n_x = 0; n_x < N128; n_x++) {
        for (int n_y = 0; n_y < N128; n_y++) {   
            fprintf(field_128_file, "%g\n", F128[n_x][n_y].value);
        }
    }
    for (int n_x = 0; n_x < N256; n_x++) {
        for (int n_y = 0; n_y < N256; n_y++) {   
            fprintf(field_256_file, "%g\n", F256[n_x][n_y].value);
        }
    }

    fclose(field_32_file);
    fclose(field_64_file);
    fclose(field_128_file);
    fclose(field_256_file);

    free(F32);
    free(F64);
    free(F128);
    free(F256);
    return 0;
}