#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void double_100(int *int_ptr) {
    for (int i = 0; i < 100; i++) {
        *(int_ptr + i) *= 2;
    }  
}


int main(void) {
    // int Array
    int array_100[100];
    for (int i = 0; i < 100; i++) {
        array_100[i] = 10 + i;
    } 
    printf("+++ Ausgabe Array +++\n");
    for (int i = 0; i < 100; i++) {
        // printf("array_100[%d] = %d\n", i, array_100[i]);
        printf("%d ", array_100[i]);
    }
    printf("\n\n");


    // int malloc
    int *int_ptr;
    int_ptr = (int*) malloc(sizeof(int) * 100);
    for (int i = 0; i < 100; i++) {
        *(int_ptr + i) = 10 + i;
    } 
    printf("+++ Ausgabe des reservierten Speicherbereichs +++\n");
    for (int i = 0; i < 100; i++) {
        // printf("*(int_ptr + [%d]) = %d\n", i, *(int_ptr + i) = 10 + i);
        printf("%d ", *(int_ptr + i));
    }
    printf("\n\n");


    // Verdoppeln der Einträge
    double_100(array_100);
    double_100(int_ptr);
    printf("+++ Ausgabe Array (verdoppelt) +++\n");
    for (int i = 0; i < 100; i++) {
        printf("%d ", array_100[i]);
    }
    printf("\n\n");
    printf("+++ Ausgabe des reservierten Speicherbereichs (verdoppelt) +++\n");
    for (int i = 0; i < 100; i++) {
        printf("%d ", *(int_ptr + i));
    }
    printf("\n\n");

    // MEME (arr[i] == *(arr + i) == *(i + arr) == i[arr])
    printf("arr[i] = %d\n", array_100[10]);
    printf("i[arr] = %d\n", 10[array_100]);

    free(int_ptr);
    return 0;
}