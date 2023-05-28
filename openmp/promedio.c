#include <stdio.h>
#include <stdlib.h>

void printMatriz(double* matriz, int N) {
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%f ", matriz[i * N + j]);
        }
        printf("\n");
    }
}

int main() {
    int* ab, * abc, * a, * p;
    int N = 4096;

    p = (int*)malloc(sizeof(int) * N * N);
    a = (int*)malloc(sizeof(int) * N * N);
    ab = (int*)malloc(sizeof(int) * N * N);
    abc = (int*)malloc(sizeof(int) * N * N);

    //inicializar matrices
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            a[i * N + j] = 0;
            p[i * N + j] = 0;
        }

    }


    //matriz a
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            a[i * N + j] += (i * N + j) + 1;
        }
    }

    //resultado ab
    int fila_I;
    for (int I = 0; I < N; I++) {
        fila_I = 0;
        for (int j = 0; j < N; j++) {
            fila_I += a[I * N + j];
        }
        for (int J = 0; J < N; J++) {
            ab[I * N + J] = fila_I;
        }
    }

    //resultado abc
    for (int I = 0; I < N; I++) {
        fila_I = 0;
        for (int j = 0; j < N; j++) {
            fila_I += ab[I * N + j];
        }
        for (int J = 0; J < N; J++) {
            abc[I * N + J] = fila_I;
        }
    }

    //resultado p
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            p[i * N + j] += abc[i * N + j] + (N * N) * abc[i * N + j];
        }
    }
    long suma = 0;
    long double promedio;
    //calcular promedio
    for (int i = 0; i < (N * N); i++) {
        suma += p[i];
    }
    promedio = (suma / (N * N));
    printf("promedio = %Lf\n", promedio);
}


