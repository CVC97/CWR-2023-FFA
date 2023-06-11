#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../cvc_numerics.h"


// Physikalische Konstanten
const double CYLINDER_HEIGHT = 0.05;                                                    // Zylinderhöhe h
const double BARRIER_HEIGHT = 0.025;                                                    // Startposition der Trennfläche y_0
const double CONCENTRATION_0 = 1e3;                                                     // Startkonzentration oberer Teil c_0_up = 1
const double DIFFUSION_CONSTANT = 4e-10;                                                // Diffusionskonstante

// Integrationsparameter
const double DELTA_Y = 1e-4;                                                            // Weite der örtlichen Diskretisierungsschritte 
const double DELTA_T = (DELTA_Y * DELTA_Y) / 2 / DIFFUSION_CONSTANT / 10;               // zeitliche Schrittweite als 1/10 des Stabilitätskriteriums
const int N = CYLINDER_HEIGHT / DELTA_Y;                                                // Anzahl der örtlichen Diskretisierungsschritte


// Diffusionsgleichung: Analytische Lösung für gegebene Position y und Zeit t
double diffusion_equation_analytical(double y, double t, void *params) {
    double D = ((double*) params)[2];                                                   // Übernahme der Diffusionskonstante aus dem params-Array
    double erf_delta_x = 1e-3;                                                          // Diskretisierung der numerischen Simpson-Integration der Fehlerfunktion
    double erf_x = (y - BARRIER_HEIGHT) / sqrt(4*D*t);                                  // Argument x der Fehlerfunktion
    return CONCENTRATION_0 / 2 * (1 - cvc_erf_simpson(erf_x, erf_delta_x));
}


// FTCS für die Iterationsvorschrift der Diffusionsgleichung
int diffusion_FTCS_step(double y[], void* params) {

    // Parameter aus dem params-Array
    double delta_t = ((double*) params)[0];
    double delta_y = ((double*) params)[1];
    double D = ((double*) params)[2];
    double h = ((double*) params)[3];

    double *y_temp = (double*) malloc(N * sizeof(double));                              // temporäres Zustands-Array für den neuen Konzentrationsschritt
    for(size_t i = 0; i < N; i++) {
        double c = y[i % N];
        double c_south = y[(i - 1) % N];
        double c_north = y[(i + 1) % N];

        if (i == 0) {                                                                   // unterer Rand: Konzentration c_i (i == 0) als Ghost
            y_temp[i] = y[i] + delta_t * D * (c - 2*c + c_north) / cvc_npow(delta_y, 2);
        } else if (i == N - 1) {                                                        // oberer Rand: Konzentration c_i (i == N - 1) als Ghost
            y_temp[i] = y[i] + delta_t * D * (c_south - 2*c + c) / cvc_npow(delta_y, 2);
        } else {                                                                        // "innere" Konzentrationen nach Iterationsvorschrift
            y_temp[i] = y[i] + delta_t * D * (c_south - 2*c + c_north) / cvc_npow(delta_y, 2);
        }
    }

    // Übertragen der neuen Werte in das Ursprungsaray
    for(int i = 0; i < N; i++) {
        y[i] = y_temp[i];
    }
    free(y_temp);
    return 0;
}


int main(void) {

    // +++++ (Aufgabe 2: Analytische Lösung der Diffusionsgleichung mit obigen Parametern) +++++

    double t_analytical[] = {cvc_HOUR, cvc_DAY, cvc_WEEK};                              // Zeitpunkte der Berechnung der analytischen Lösung (1 Stunde, 1 Tag, 1 Woche)
    double parameters[] = {DELTA_T, DELTA_Y, DIFFUSION_CONSTANT, CYLINDER_HEIGHT};      // Array mit Parametern zur Übergabe
    double y, t, c_y_t;                                                                 // Variablen der Position y und Zeit t, sowie analytischer Wert der Konzentration

    // Datei der analytischen Werte mit y-Positionen als Zeile 1
    FILE* diffusion_analyt_file = fopen("A23_Diffusion_Analytical.csv", "w");                
    fprintf(diffusion_analyt_file, "0.05");
    for (int i_y = 0; i_y < N; i_y++) {
        y = (i_y+1) * DELTA_Y;
        fprintf(diffusion_analyt_file, ", %g", 0.05-y);
    }

    // Beschreiben der Datei mit den analytischen Werten der 3 Zeiten
    for (int i_t = 0; i_t < 3; i_t++) {                                                 // Itieren über die 3 Zeiten
        t = t_analytical[i_t];
        c_y_t = diffusion_equation_analytical(5, t, parameters);
        fprintf(diffusion_analyt_file, "\n%g", c_y_t);
        for (int i_y = 0; i_y < N; i_y++) {                                             // Itieren über die N y-Positionen
            y = (i_y+1) * DELTA_Y;
            c_y_t = diffusion_equation_analytical(y, t, parameters);
            fprintf(diffusion_analyt_file, ", %g", c_y_t);
        }
    }
    fclose(diffusion_analyt_file);


    // +++++ (Aufgabe 5: Numerische Lösung der Diffusionsgleichung für die 3 bekannten Zeiten) +++++

    // Speicher der Zustandarrays für die 3 Zeiten
    double *y_hour = (double*) malloc((N+1) * sizeof(double));
    double *y_day = (double*) malloc((N+1) * sizeof(double)); 
    double *y_week = (double*) malloc((N+1) * sizeof(double));  

    // Dateien der numerischen Werte
    FILE* diffusion_numeric_hour_file = fopen("A23_Diffusion_Numerical_Hour.csv", "w");
    FILE* diffusion_numeric_day_file = fopen("A23_Diffusion_Numerical_Day.csv", "w");
    FILE* diffusion_numeric_week_file = fopen("A23_Diffusion_Numerical_Week.csv", "w");

    // Initialisieren der Zustandsarrays mit den Startwerten
    for (int i_y = 0; i_y < N+1; i_y++) {
        y = i_y * DELTA_Y;
        if (y <= BARRIER_HEIGHT) {                                                      // Unter / "im" dem Plättchens: Konzentration gleich 0
            y_hour[i_y] = 0;
            y_day[i_y] = 0;
            y_week[i_y] = 0;
        } else {                                                                        // Über dem Plättchen: Ausgangskonzentration
            y_hour[i_y] = CONCENTRATION_0;
            y_day[i_y] = CONCENTRATION_0;
            y_week[i_y] = CONCENTRATION_0;
        }
    }

    // Numerische Integration bis t gleich 1 Stunde
    t = 0;
    while (t < cvc_HOUR) {
        diffusion_FTCS_step(y_hour, parameters);
        diffusion_FTCS_step(y_day, parameters);
        diffusion_FTCS_step(y_week, parameters);
        t += DELTA_T;
    }
    for (int i_y = 0; i_y < N+1; i_y++) {                                               // Zustand des Systems nach 1 Stunde                                                                             
        fprintf(diffusion_numeric_hour_file, "%g\n", y_hour[i_y]);
    }
    fclose(diffusion_numeric_hour_file);

    // Numerische Integration bis t gleich 1 Tag
    while (t < cvc_DAY) {
        diffusion_FTCS_step(y_day, parameters);
        diffusion_FTCS_step(y_week, parameters);
        t += DELTA_T;
    }
    for (int i_y = 0; i_y < N+1; i_y++) {                                               // Zustand des Systems nach 1 Tag       
        fprintf(diffusion_numeric_day_file, "%g\n", y_day[i_y]);
    }
    fclose(diffusion_numeric_day_file);

    // Numerische Integration bis t gleich 1 Woche
    while (t < cvc_WEEK) {
        diffusion_FTCS_step(y_week, parameters);
        t += DELTA_T;
    }
    for (int i_y = 0; i_y < N+1; i_y++) {                                               // Zustand des Systems nach 1 Woche       
        fprintf(diffusion_numeric_week_file, "%g\n", y_week[i_y]);
    }
    fclose(diffusion_numeric_week_file);
    return 0;
}