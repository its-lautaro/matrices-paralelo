/* Producto punto de matrices
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
int BS, bloques_por_thread;
double* M1, * M2, * Res;

void printMatriz(double* matriz, int n) {
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%f ", matriz[i * n + j]);
        }
        printf("\n");
    }
}

void* dot_prod(void* arg) {
    int id = *(int*)arg;

    for (int bloque = 0; bloque < bloques_por_thread; bloque++) {
        for (int fila = BS * bloque;fila < (BS * bloque) + BS;fila++) {
            for (int columna = id * BS; columna < (id * BS) + BS; columna++) {
                //consultar!
                //Res[fila * N + columna] += M1[fila * N + columna] * M2[fila * N + columna];
                //printf("%d: Bloque %d Matriz resultado[%d,%d]=M1[%d,%d]+M2[%d,%d]\n",id, bloque, fila, columna, fila, columna, columna, fila);
            }
        }
    }

    pthread_exit(NULL);
}
 
void dot_prod_paralelo() {
    pthread_t threads[T];
    int ids[T];

    for (int i = 0; i < T; i++) {
        ids[i] = i;
        pthread_create(&threads[i], NULL, &dot_prod, (void*)&ids[i]);
    }
    for (int i = 0; i < T; i++) {
        pthread_join(threads[i], NULL);
    }
}

int main(int argc, char const* argv[]) {
    N = 8;
    BS = 2;
    T = 4;
    M1 = (double*)malloc(sizeof(double) * N * N);
    M2 = (double*)malloc(sizeof(double) * N * N);
    Res = (double*)malloc(sizeof(double) * N * N);

    bloques_por_thread = ((N * N) / (BS * BS)) / T;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            M1[i * N + j] = 1.0;
            M2[i * N + j] = 1.0;
            Res[i * N + j] = 0.0;
        }
    }

    dot_prod_paralelo();

    printMatriz(Res, N);

    return 0;
}