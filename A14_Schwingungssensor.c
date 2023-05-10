#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cvc_numerics.h"


// Physikalische Konstanten des Systems
const int N = 2;                                        // Dimension des Problems N
const double mass = 1;                                  // Masse mass
const double k = 10;                                    // Federkonstante k beider Federn           
const double r_0 = 1;                                   // Ruhelänge r_0 beider Federn
const double r_A[2] = {-1, 0};                          // Position der Verankerung der linken Feder r_A
const double r_B[2] = {1, 0};                           // Position der Verankerung der rechten Feder r_b

const double smoothing_factor = 1e-6;                   // Vermeidung von numerischen Fehlern bei Division (durch 0)


// (gewöhnliche) Differentialgleichung des Federsystems
int ODE_dual_springs(double t, const double y[], double f[], void *params) {
    // Übertragen der Geschwindigkeiten in das Ableitungsarray
    f[0] = y[2];
    f[1] = y[3];

    // Berechnung der Abstände der Masse zur jeweiligen Federverankerung
    double dist_r_A = norm_2D(y[0] - r_A[0], y[1] - r_A[1]) + smoothing_factor;
    double dist_r_B = norm_2D(y[0] - r_B[0], y[1] - r_B[1]) + smoothing_factor;

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
    int dimension = 2 * N;
    double T_x = 2 * CVC_PI / sqrt(2*k);
    double x_0 = -0.5;
    double t = 0, delta_t = 1e-4;


// +++++ (Aufgabe 4, 5: Position x(t) in Abhängigkeit der Zeit) +++++

    // Initialisierung der Startposition der Masse
    double y_euler[4] = {-0.5, 0}; 
    double y_rk4[4] = {-0.5, 0};

    // Ausgabedatei
    FILE* pos_file = fopen("data/A14_Schwingungssensor_pos.csv", "w"); 
    fprintf(pos_file, "Zeit, x(t) [Euler], y(t) [Euler], x(t) [RK4], y(t) [RK4]\n");
    fprintf(pos_file, "%g, %g, %g, %g, %g\n", t, y_euler[0], y_euler[1], y_rk4[0], y_rk4[1]);

     // Integration bis 10 * T_x
    while (t <= 10 * T_x) {
        t += delta_t;
        euler_step(t, delta_t, y_euler, ODE_dual_springs, dimension, NULL);
        rk4_step(t, delta_t, y_rk4, ODE_dual_springs, dimension, NULL);
        fprintf(pos_file, "%g, %g, %g, %g, %g\n", t, y_euler[0], y_euler[1], y_rk4[0], y_rk4[1]);
    }
    fclose(pos_file);


// +++++ (Aufgabe 6: Abweichung der numerischen von der analytischen Lösung für variierte Zeitschritte delta_t) +++++

    // Ausgabdatei
    FILE* res_file = fopen("data/A14_Schwingungssensor_res.csv", "w");
    fprintf(res_file, "delta_t, Residue Euler, Residue RK4\n"); 

    // Itereation über 100 logarithmisch gleichverteilte delta_t
    for (int i = 0; i < 100; i++) {
        t = 0, delta_t = pow(10, -1-4*i/99.0);          // Neuberechnung der Zeitparameter für jede Integration

        // (Re)initialisierung der Startposition der Masse
        double y_euler_res[4] = {-0.5, 0}; 
        double y_rk4_res[4] = {-0.5, 0};

        // Integration bis 10 * T_x
        while (t + delta_t < 10 * T_x) {
            t += delta_t;
            euler_step(t, delta_t, y_euler_res, ODE_dual_springs, dimension, NULL);
            rk4_step(t, delta_t, y_rk4_res, ODE_dual_springs, dimension, NULL);
        }
        
        // Verwendung eines zusätzlichen Zeitschritts (Differenz zur analytischen T_x) für Vergleichbarkeit mit analytischer Lösung
        euler_step(t, 10*T_x - t, y_euler_res, ODE_dual_springs, dimension, NULL);
        rk4_step(t, 10*T_x - t, y_rk4_res, ODE_dual_springs, dimension, NULL);

        // Ausgeben der Residuen
        fprintf(res_file, "%g, %g, %g\n", delta_t, fabs(y_euler_res[0] - x_0), fabs(y_rk4_res[0] - x_0));
    }
    fclose(res_file);


// +++++ (Aufgabe 7: Numerische Berechnung der Trajektorie mit RK4 für drei unterschiedliche Ausgangszustände) +++++
    t = 0, delta_t = 10e-4;

    // Initialisierung der Startpositionen der Masse
    double y1_rk4[4] = {-0.5, 0.1, 0.4};
    double y2_rk4[4] = {-0.5, 0.2, 0, -0.1};
    double y3_rk4[4] = {0, -0.3, 1};

    // Ausgabedatei
    FILE* fancy_pos_file = fopen("data/A14_Schwingungssensor_fancy_pos.csv", "w"); 
    fprintf(fancy_pos_file, "Zeit, x1(t), y1(t), x2(t), y2(t), x3(t), y3(t)\n");
    fprintf(fancy_pos_file, "%g, %g, %g, %g, %g, %g, %g\n", t, y1_rk4[0], y1_rk4[1], y2_rk4[0], y2_rk4[1], y3_rk4[0], y3_rk4[1]);

     // Integration bis 10 * T_x
    while (t <= 10 * T_x) {
        t += delta_t;
        rk4_step(t, delta_t, y1_rk4, ODE_dual_springs, dimension, NULL);
        rk4_step(t, delta_t, y2_rk4, ODE_dual_springs, dimension, NULL);
        rk4_step(t, delta_t, y3_rk4, ODE_dual_springs, dimension, NULL);
        fprintf(fancy_pos_file, "%g, %g, %g, %g, %g, %g, %g\n", t, y1_rk4[0], y1_rk4[1], y2_rk4[0], y2_rk4[1], y3_rk4[0], y3_rk4[1]);
    }
    fclose(fancy_pos_file);
    return 0;
}