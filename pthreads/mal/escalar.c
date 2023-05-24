/* Producto escalar de matrices
* Memoria Compartida (pthread)
* Lautaro
*/
#include<stdio.h>
#include<stdlib.h>  
#include <sys/time.h>
#include <pthread.h>

//variables compartidas
int N;//tamaño de la matriz
int T;//cantidad de threads
int BS;
double* matriz;
double escalar;
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

void* producto_escalar(void* args) {
    int id = *((int*)args);
    for (int bloque = 0; bloque < bloques_por_thread; bloque++) {
        for (int fila = BS * bloque;fila < (BS * bloque) + BS;fila++) {
            for (int columna = id * BS; columna < (id * BS) + BS; columna++) {
                matriz[fila * N + columna] *= escalar;
            }
        }
    }


    pthread_exit(NULL);
}

void producto_escalar_paralelo() {
    pthread_t threads[T];
    int ids[T];

    for (int i = 0; i < T; i++) {
        ids[i] = i;
        pthread_create(&threads[i], NULL, &producto_escalar, (void*)&ids[i]);
    }
    for (int i = 0; i < T; i++) {
        pthread_join(threads[i], NULL);
    }
}

int main(int argc, char const* argv[]) {
    N = 8;
    BS = 2;
    T = 4;
    matriz = (double*)malloc(sizeof(double) * N * N);
    escalar = 3.0;

    bloques_por_thread = ((N * N) / (BS * BS)) / T;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matriz[i * N + j] = 1.0;
        }
    }

    producto_escalar_paralelo();

    printMatriz(matriz, N);

    return 0;
}