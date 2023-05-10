#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// Struktur eines 2er Tupels
struct tupel_2 {
    double sol1;
    double sol2;
};


// Analytische Lösungen der gegebenen Gleichung
const double sol1_analytical = 1e8;
const double sol2_analytical = 1e-8;


// Lösung der Mitternachts-Formel
struct tupel_2 MiNaFo(double a, double b, double c) {
    double sol1, sol2;

    sol1 = (-b + sqrt(pow(b, 2) - 4*a*c)) / (2*a);
    sol2 = (-b - sqrt(pow(b, 2) - 4*a*c)) / (2*a);

    struct tupel_2 minafo;
    minafo.sol1 = sol1;
    minafo.sol2 = sol2;
    return minafo;
}


// Lösung der erweiterten Mitternachts-Formel
struct tupel_2 MiNaFo_extended(double a, double b, double c) {
    double sol1, sol2;
    sol1 = (2*c) / (-b - sqrt(pow(b, 2) - 4*a*c));
    sol2 = (2*c) / (-b + sqrt(pow(b, 2) - 4*a*c));

    struct tupel_2 minafo_extended;
    minafo_extended.sol1 = sol1;
    minafo_extended.sol2 = sol2;
    return minafo_extended;
}


// kombinierte Lösungsmethode quadratischer Gleichungen
struct tupel_2 solve_quadratic(double a, double b, double c) {
    double sol1, sol2;

    if (b > 0) {
        sol1 = (2*c) / (-b - sqrt(pow(b, 2) - 4*a*c));
        sol2 = (-b - sqrt(pow(b, 2) - 4*a*c)) / (2*a);
    } else {
        sol1 = (-b + sqrt(pow(b, 2) - 4*a*c)) / (2*a);
        sol2 = (2*c) / (-b + sqrt(pow(b, 2) - 4*a*c));
    }

    struct tupel_2 quadratic_solution;
    quadratic_solution.sol1 = sol1;
    quadratic_solution.sol2 = sol2;
    return quadratic_solution;
}


int main(void) {
    double a, b, c;
    a = 1;
    b = - (1e16 + 1) / 1e8;
    c =  1;
    
    struct tupel_2 minafo, minafo_extended, quadratic_solution;
    minafo = MiNaFo(a, b, c);
    minafo_extended = MiNaFo_extended(a, b, c);
    quadratic_solution = solve_quadratic(a, b, c);

    double minafo_rel_error1, minafo_rel_error2, minafo_extended_rel_error1, minafo_extended_rel_error2, quadratic_solution_rel_error1, quadratic_solution_rel_error2;
    minafo_rel_error1 = (sol1_analytical - minafo.sol1) / sol1_analytical * 100;
    minafo_rel_error2 = (sol2_analytical - minafo.sol2) / sol2_analytical * 100;

    minafo_extended_rel_error1 = (sol1_analytical - minafo_extended.sol1) / sol1_analytical * 100;
    minafo_extended_rel_error2 = (sol2_analytical - minafo_extended.sol2) / sol2_analytical * 100;

    quadratic_solution_rel_error1 = (sol1_analytical - quadratic_solution.sol1) / sol1_analytical * 100;
    quadratic_solution_rel_error2 = (sol2_analytical - quadratic_solution.sol2) / sol2_analytical * 100;    

    printf("+++ MiNaFo +++\nx1 = %g (Relativer Fehler: %g %%)\nx2 = %g (%g %%)\n", minafo.sol1, minafo_rel_error1, minafo.sol2, minafo_rel_error2);
    printf("\n+++ Erweiterne MiNaFo +++\nx1 = %g (%g %%)\nx2 = %g (%g %%)\n", minafo_extended.sol1, minafo_extended_rel_error1, minafo_extended.sol2, minafo_extended_rel_error2);
    printf("\n+++ KING (I have no idea what I'm doing) +++\nx1 = %g (%g %%)\nx2 = %g (%g %%)\n", quadratic_solution.sol1, quadratic_solution_rel_error1, quadratic_solution.sol2, quadratic_solution_rel_error2);
    return 0;
}