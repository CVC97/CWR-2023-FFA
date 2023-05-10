#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// Populationskonstanten
const double n_0 = 0.5;           // Anfangsgröße der Population
const double r_0 = 4;             // Reproduktionsrate
const double K = 4;               // Kapazität des Lebensraumes


// Berechnet für gegebenes my und x(i) den nächsten Schritt x(i+1)
double x_i1(double mu, double x_i) {
    return 4*mu*x_i*(1-x_i);
}


// Berechnet die analytische Lösung
double n(double t) {
    return n_0*exp(r_0*t) / (1 + K*n_0 * (exp(r_0*t - 1)));
}


int main(void) {
    // Stabilitätsanalyse Population(t)
    double mu_1 = 0.4;
    double mu_2 = 0.74;
    double mu_3 = 0.77;

    // Ausgangswerte 
    double x_mu_1 = n_0;
    double x_mu_2 = n_0;
    double x_mu_3 = n_0;
    double x_analyt = n_0;

    double t_1 = 0;
    double t_2 = 0;
    double t_3 = 0;

    // Ausgabe-Dateien
    FILE* x_file = fopen("data/A11_Populationsdynamik_t.csv", "w");
    fprintf(x_file, "Zeit t_1, x_1, t_2, x_2, t_3, x_3, n_analytisch\n");
    fprintf(x_file, "%g, %g, %g, %g, %g, %g, %g\n", t_1, t_2, t_3, x_mu_1, x_mu_2, x_mu_3, x_analyt);

    // Integrationsschleife x(t)
    for (int i = 0; i < 100; i++) {
        t_1 += mu_1;
        t_2 += mu_2;
        t_3 += mu_3;

        x_mu_1 = x_i1(mu_1, x_mu_1);
        x_mu_2 = x_i1(mu_2, x_mu_2);
        x_mu_3 = x_i1(mu_3, x_mu_3);
        x_analyt = n(t_3);

        fprintf(x_file, "%g, %g, %g, %g, %g, %g, %g\n", t_1, t_2, t_3, x_mu_1, x_mu_2, x_mu_3, x_analyt);
    }
    fclose(x_file);


    // Stabilitätsanalyse Populationsgröße nach 1000 Integrationsschritten abhängig von my
    double mu = 0.4;

    // Ausgabe-Dateiein x_1000(mu)
    FILE* x1000_file = fopen("data/A11_Populationsdynamik_mu.csv", "w");
    fprintf(x1000_file, "mu, x_1000\n");
    while (mu < 1) {
        double x_1000 = n_0;
        for (int i = 0; i < 1000; i++) {
            x_1000 = x_i1(mu, x_1000);
        }
        fprintf(x1000_file, "%g, %g\n", mu, x_1000);
        mu += 0.01;
    }
    fclose(x1000_file);
    return 0;
}