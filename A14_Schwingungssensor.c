#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../cvc_numerics.h"


// Physikalische Konstanten des Systems
const int N = 2;                                                                // Dimension N des Problems
const double mass = 1;                                                          // Masse mass
const double k = 10;                                                            // Federkonstante k beider Federn           
const double r_0 = 1;                                                           // Ruhelänge r_0 beider Federn
const double r_A[2] = {-1, 0};                                                  // Position der Verankerung der linken Feder r_A
const double r_B[2] = {1, 0};                                                   // Position der Verankerung der rechten Feder r_b

const double smoothing_factor = 1e-12;                                          // Vermeidung von numerischen Fehlern bei Division (durch 0)


// (gewöhnliche) Differentialgleichung des Federsystems
int ODE_dual_springs(double t, const double y[], double f[], void *params) {
    // Übertragen der Geschwindigkeiten in das Ableitungsarray
    f[0] = y[2];
    f[1] = y[3];

    // Berechnung der Abstände der Masse zur jeweiligen Federverankerung
    double dist_r_A = cvc_norm_2D(y[0] - r_A[0], y[1] - r_A[1]) + smoothing_factor;
    double dist_r_B = cvc_norm_2D(y[0] - r_B[0], y[1] - r_B[1]) + smoothing_factor;

    // Berechnung der richtungslosen Kräfte beider Federn
    double f_A = -k * (dist_r_A - r_0);
    double f_B = -k * (dist_r_B - r_0);

    // Berechnung der x- und y-Komponenten der Kräfte beider Federn
    double f_A_x = f_A * (y[0] - r_A[0]) / dist_r_A;
    double f_A_y = f_A * (y[1] - r_A[1]) / dist_r_A;

    double f_B_x = f_B * (y[0] - r_B[0]) / dist_r_B;
    double f_B_y = f_B * (y[1] - r_B[1]) / dist_r_B;

    // Übertragen / Berechnen der x- und y-Komponenten der insgesamt wirkenden Kraft
    f[2] = (f_A_x + f_B_x) / mass;
    f[3] = (f_A_y + f_B_y) / mass;
    return 0;
}


int main(void) {
    // Simulationsparameter
    int dimension = 2 * N;                                                      // Dimensionalität des Zustandsvektors des Systems
    double T_x = 2 * cvc_PI / sqrt(2*k), T_x_Euler = 0, T_x_RK4 = 0;            // Periodendauern T_x
    double x_0 = -0.5;                                                          // Startkoordinate x_0
    double t = 0, delta_t = 1e-4;                                               // Startzeit t = 0, Größe der Zeitschritte delta_t


// +++++ (Aufgabe 4, 5: Position x(t) in Abhängigkeit der Zeit) +++++

    // Initialisierung der Startposition der Masse
    double y_euler[4] = {x_0, 0}; 
    double y_rk4[4] = {x_0, 0};

    // numerische Nullstellensuche über Vorzeichenwechsel der x-Koordinate
    int euler_number_zero_crossings = 0;                                        // Anzahl der Nulldurchgänge (Euler)
    int rk4_number_zero_crossings = 0;                                          // Anzahl der Nulldurchgänge (RK4)
    double *euler_zero_crossings = (double*) malloc(sizeof(double));            // Reservierung Speicher für Nulldurchgänge (Euler)
    double *rk4_zero_crossings = (double*) malloc(sizeof(double));              // Reservierung Speicher für Nulldurchgänge (RK4)

    // Ausgabedatei
    FILE* pos_file = fopen("A14_Schwingungssensor_pos.csv", "w"); 
    fprintf(pos_file, "Zeit, x(t) [Euler], y(t) [Euler], x(t) [RK4], y(t) [RK4]\n");
    fprintf(pos_file, "%g, %g, %g, %g, %g\n", t, y_euler[0], y_euler[1], y_rk4[0], y_rk4[1]);

    // Integration über 10 Periodendauern T_x
    while (t <= 10 * T_x) {
        t += delta_t;

        double x_pre_euler = y_euler[0];                                        // x-Koordinate vor Integration (Euler) 
        double x_pre_rk4 = y_rk4[0];                                            // x-Koordinate vor Integration (RK4)

        cvc_euler_step(t, delta_t, y_euler, ODE_dual_springs, dimension, NULL);
        cvc_rk4_step(t, delta_t, y_rk4, ODE_dual_springs, dimension, NULL);

        double x_post_euler = y_euler[0];                                       // x-Koordinate nach Integration (Euler) 
        double x_post_rk4 = y_rk4[0];                                           // x-Koordinate nach Integration (Euler) 

        // Bestimmung der Zeiten der Nulldurchgänge
        if (x_pre_euler * x_post_euler < 0) {                                   // Prüfen eines Nulldurchganges über Vorzeichenwechsel (Euler)         
            euler_number_zero_crossings++;               
            euler_zero_crossings = (double*) realloc(euler_zero_crossings, sizeof(double) * euler_number_zero_crossings);
            euler_zero_crossings[euler_number_zero_crossings - 1] = t;          // Eintragen des Nulldurchganges in vergrößerten Speicherbereich (Euler)             

        }
        if (x_pre_rk4 * x_post_rk4 < 0) {                                       // Prüfen eines Nulldurchganges über Vorzeichenwechsel (RK4)
            rk4_number_zero_crossings++;
            rk4_zero_crossings = (double*) realloc(rk4_zero_crossings, sizeof(double) * rk4_number_zero_crossings);
            rk4_zero_crossings[rk4_number_zero_crossings - 1] = t;              // Eintragen des Nulldurchganges in vergrößerten Speicherbereich (RK4) 

        }
        fprintf(pos_file, "%g, %g, %g, %g, %g\n", t, y_euler[0], y_euler[1], y_rk4[0], y_rk4[1]);
    }

    // Ausgabe der Periodendauer (wenn möglich) 
    printf("Periodendauer T_x (analytisch): %.14g\n", T_x);
    if (euler_number_zero_crossings < 2) {                                      // Prüfen ob ausreichend Nulldurchgänge (Euler)
        printf("ERROR! Insufficient number of zero-crossings for Euler.");
    } else {                                                                    // Berechnung T_x für ausreichend Nulldurchgänge (Euler)
        T_x_Euler = 2.0 / (euler_number_zero_crossings - 1) * (euler_zero_crossings[euler_number_zero_crossings - 1] - euler_zero_crossings[0]);
        printf("Periodendauer T_x (Euler): %.14g\n", T_x_Euler);                // Printen Periodendauer als 2 * [t(Nulldurchgang N) - t(Nulldurchgang 1)] / (N - 1) für Nulldurchgänge > 1 (Euler)
    }
    if (rk4_number_zero_crossings < 2) {                                        // Prüfen ob ausreichend Nulldurchgänge (RK4)
        printf("ERROR! Insufficient number of zero-crossings for RK4.");
        return 1;
    } else {                                                                    // Berechnung T_x für ausreichend Nulldurchgänge (RK4)                           
        T_x_RK4 = 2.0 / (rk4_number_zero_crossings - 1) * (rk4_zero_crossings[rk4_number_zero_crossings - 1] - rk4_zero_crossings[0]);
        printf("Periodendauer T_x (RK4): %.14g\n", T_x_RK4);                    // Printen Periodendauer als 2 * [t(Nulldurchgang N) - t(Nulldurchgang 1)] / (N - 1) für Nulldurchgänge > 1 (RK4)  
    }
    fclose(pos_file), free(euler_zero_crossings), free(rk4_zero_crossings);


// +++++ (Aufgabe 6: Abweichung der numerischen von der analytischen Lösung für variierte Zeitschritte delta_t) +++++

    // Ausgabdatei
    FILE* res_file = fopen("A14_Schwingungssensor_res.csv", "w");
    fprintf(res_file, "delta_t, Residue Euler, Residue RK4\n"); 

    // Itereation über 100 logarithmisch gleichverteilte Zeitschritte delta_t
    for (int i = 0; i < 100; i++) {
        t = 0, delta_t = pow(10, -1-4*i/99.0);                                  // Neuberechnung der Zeitparameter für jede Integration

        // (Re)initialisierung der Startposition der Masse
        double y_euler_res[4] = {x_0, 0}; 
        double y_rk4_res[4] = {x_0, 0};

        // Integration über 10 Periodendauern T_x_RK4
        while (t + delta_t < 10 * T_x_RK4) {
            t += delta_t;
            cvc_euler_step(t, delta_t, y_euler_res, ODE_dual_springs, dimension, NULL);
            cvc_rk4_step(t, delta_t, y_rk4_res, ODE_dual_springs, dimension, NULL);
        }
        
        // zusätzlicher Zeitschritt (Differenz zur analytischen T_x_RK4) für Vergleichbarkeit mit analytischer Lösung
        cvc_euler_step(t, 10*T_x_RK4 - t, y_euler_res, ODE_dual_springs, dimension, NULL);
        cvc_rk4_step(t, 10*T_x_RK4 - t, y_rk4_res, ODE_dual_springs, dimension, NULL);

        // Ausgeben der Residuen
        fprintf(res_file, "%g, %g, %g\n", delta_t, fabs(y_euler_res[0] - x_0), fabs(y_rk4_res[0] - x_0));
    }
    fclose(res_file);


// +++++ (Aufgabe 7: Numerische Berechnung der Trajektorie mit RK4 für drei unterschiedliche Ausgangszustände) +++++
    t = 0, delta_t = 10e-4;

    // Initialisierung der Startpositionen der Masse
    double y1_rk4[4] = {-0.75, 0.2, 0, -0.2};
    double y2_rk4[4] = {-0.5, -0.5};
    double y3_rk4[4] = {0, -0.3, 1, 1};

    // Ausgabedatei
    FILE* fancy_pos_file = fopen("A14_Schwingungssensor_fancy_pos.csv", "w"); 
    fprintf(fancy_pos_file, "Zeit, x1(t), y1(t), x2(t), y2(t), x3(t), y3(t)\n");
    fprintf(fancy_pos_file, "%g, %g, %g, %g, %g, %g, %g\n", t, y1_rk4[0], y1_rk4[1], y2_rk4[0], y2_rk4[1], y3_rk4[0], y3_rk4[1]);

    // Integration über 10 Periodendauern T_x_RK4
    while (t <= 10 * T_x_RK4) {
        t += delta_t;
        cvc_rk4_step(t, delta_t, y1_rk4, ODE_dual_springs, dimension, NULL);
        cvc_rk4_step(t, delta_t, y2_rk4, ODE_dual_springs, dimension, NULL);
        cvc_rk4_step(t, delta_t, y3_rk4, ODE_dual_springs, dimension, NULL);
        fprintf(fancy_pos_file, "%g, %g, %g, %g, %g, %g, %g\n", t, y1_rk4[0], y1_rk4[1], y2_rk4[0], y2_rk4[1], y3_rk4[0], y3_rk4[1]);
    }
    fclose(fancy_pos_file);
    return 0;
}