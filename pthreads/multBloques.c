/* Realiza los productos A*B*C, D*C*B y calcula maximo y minimo
* Memoria Compartida (pthread)
* Lautaro
*/
#include<stdio.h>
#include<stdlib.h>  
#include <sys/time.h>
#include <pthread.h>

//variables compartidas
int N;//tamaño de la matriz
int BS;//tamaño de bloque
int T;//cantidad de threads
double* A, * B, * AB, * C, * ABC, * D, * DC, * DCB;
double maxD = -1, minA = 9999;

void printMatriz(double* matriz, int n) {
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%f ", matriz[i * n + j]);
        }
        printf("\n");
    }
}

double dwalltime() {
    double sec;
    struct timeval tv;

    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

void* multBloques(void* arg) {
    int id = *(int*)arg;
    double* Ablk, * Bblk, * Cblk, * Dblk, * ABblk, * ABCblk, * DCblk, * DCBblk;
    for (int I = (id * BS); I < N; I += (T * BS)) {
        for (int J = 0; J < N; J += BS) {
            ABblk = &AB[I * N + J];
            DCblk = &DC[I * N + J];
            for (int K = 0; K < N; K += BS) {
                Ablk = &A[I * N + K];
                Bblk = &B[J * N + K];
                Dblk = &D[I * N + K];
                Cblk = &C[J * N + K];
                for (int i = 0; i < BS; i++) {
                    for (int j = 0; j < BS; j++) {
                        for (int k = 0; k < BS; k++) {
                            //AB=A*B
                            ABblk[i * N + j] += Ablk[i * N + k] * Bblk[j * N + k];
                            //DC=D*C
                            DCblk[i * N + j] += Dblk[i * N + k] * Cblk[j * N + k];

                        }
                        //MIN(A)
                        if (Ablk[i * N + j] < minA) minA = Ablk[i * N + j];
                        //MAX(D)
                        if (Dblk[i * N + j] > maxD) maxD = Dblk[i * N + j];
                    }
                }
            }
        }
        for (int J = 0; J < N; J += BS) {
            ABCblk = &ABC[I * N + J];
            DCBblk = &DCB[I * N + J];
            for (int K = 0; K < N; K += BS) {
                ABblk = &AB[I * N + K];
                Cblk = &C[J * N + K];
                DCblk = &DC[I * N + K];
                Bblk = &B[J * N + K];
                for (int i = 0; i < BS; i++) {
                    for (int j = 0; j < BS; j++) {
                        for (int k = 0; k < BS; k++) {
                            //ABC=AB*C
                            ABCblk[i * N + j] += ABblk[i * N + k] * Cblk[j * N + k];
                            //DCB=DC*B
                            DCBblk[i * N + j] += DCblk[i * N + k] * Bblk[j * N + k];
                        }
                    }
                }
            }
        }
    }

    pthread_exit(NULL);
}

int main(int argc, char* argv[]) {
    BS = 512;
    N = 1024;
    T = 4;

    double timetick;

    pthread_t hilos[T];
    int threads_ids[T];

    // Alocar  
    A = (double*)malloc(N * N * sizeof(double));
    B = (double*)malloc(N * N * sizeof(double));
    AB = (double*)malloc(N * N * sizeof(double));
    C = (double*)malloc(N * N * sizeof(double));
    ABC = (double*)malloc(N * N * sizeof(double));
    D = (double*)malloc(N * N * sizeof(double));
    DC = (double*)malloc(N * N * sizeof(double));
    DCB = (double*)malloc(N * N * sizeof(double));

    // Inicializacion
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i * N + j] = 1;
            B[j * N + i] = 1.0; //ordenada por columnas
            AB[i * N + j] = 0.0;
            C[j * N + i] = 1.0; //ordenada por columnas
            ABC[i * N + j] = 0.0;
            D[i * N + j] = 1;
            DC[i * N + j] = 0.0;
            DCB[i * N + j] = 0.0;
        }
    }


    //tomar tiempo start
    timetick = dwalltime();

    for (int i = 0; i < T; i++) {
        threads_ids[i] = i;
        pthread_create(&hilos[i], NULL, &multBloques, (void*)&threads_ids[i]);
    }
    for (int i = 0; i < T; i++) {
        pthread_join(hilos[i], NULL);
    }

    double totalTime = dwalltime() - timetick;
    printf("Tiempo: %fs\n", totalTime);

    //estoy multiplicando bien?
    // printf("Matriz AB\n");
    // printMatriz(AB, N);
    // printf("Matriz DC\n");
    // printMatriz(DC, N);
    // printf("Matriz ABC\n");
    // printMatriz(ABC, N);
    // printf("Matriz DCB\n");
    // printMatriz(DCB, N);

    // printf("MaximoD: %f, MinimoA:%f\n", maxD, minA);
    return 0;
}