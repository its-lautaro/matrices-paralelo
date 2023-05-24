/* Suma por bloques
* Memoria Compartida (pthread)
*
*/
#include<stdio.h>
#include<stdlib.h>  
#include <sys/time.h>
#include <pthread.h>

//variables compartidas
int N;//tamaño de la matriz
int T;//cantidad de threads
int BS;
double* ABC, * DCB, * P;
int bloques_por_thread;

void printMatriz(double* matriz, int n) {
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%f ", matriz[i * n + j]);
        }
        printf("\n");
    }
}

void* suma(void* args) {
    int id = *((int*)args);
    for (int bloque = 0; bloque < bloques_por_thread; bloque++) {
        for (int fila = BS * bloque;fila < (BS * bloque) + BS;fila++) {
            for (int columna = id * BS; columna < (id * BS) + BS; columna++) {
                //             matriz[fila * N + columna] *= escalar;
                P[fila*N+columna] = ABC[fila*N+columna]+DCB[fila*N+columna];
            }
        }
    }


    pthread_exit(NULL);
}

void suma_paralelo() {
    pthread_t threads[T];
    int ids[T];

    for (int i = 0; i < T; i++) {
        ids[i] = i;
        pthread_create(&threads[i], NULL, &suma, (void*)&ids[i]);
    }
    for (int i = 0; i < T; i++) {
        pthread_join(threads[i], NULL);
    }
}

int main(int argc, char const* argv[]) {
    N = 8;
    BS = 2;
    T = 4;
    ABC = (double*)malloc(sizeof(double) * N * N);
    DCB = (double*)malloc(sizeof(double) * N * N);
    P = (double*)malloc(sizeof(double) * N * N);

    bloques_por_thread = ((N * N) / (BS * BS)) / T;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            ABC[i * N + j] = 32.0;
            DCB[i * N + j] = 58.0;
        }
    }

    suma_paralelo();

    printMatriz(P, N);

    return 0;
}