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
double* maximos, * minimos;
double maxD, minA;

pthread_barrier_t barrera;

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
    double local_min = 9999, local_max = -1;
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
                        if (Ablk[i * N + j] < local_min) {
                            local_min = Ablk[i * N + j];
                            //printf("proceso %d local_min %f\n", id, local_min);
                        }

                        //MAX(D)
                        if (Dblk[i * N + j] > local_max) {
                            local_max = Dblk[i * N + j];
                            //printf("proceso %d local_max %f\n", id, local_max);
                        }
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
    //cargo minimo y maximo local
    minimos[id] = local_min;
    maximos[id] = local_max;
    pthread_barrier_wait(&barrera); //esperar a que se terminen de escribir las variables compartidas
    //reducir arreglo
    int offset = 2;
    while (offset <= T) {
        if (id % offset == 0) {
            local_min = 9999;
            local_max = -1;
            for (int i = (id/offset);i < ((id)/offset)+offset;i++) { //cada hilo recorre su parte del arreglo y busca los minimos y maximos
                if (minimos[i] < local_min) local_min = minimos[i];
                if (maximos[i] > local_max) local_max = maximos[i];
            }
        }
        pthread_barrier_wait(&barrera); //esperan todos a que se termine de leer las variables compartidas
        if (id % offset == 0) {
            //se reduce el arreglo
            maximos[id / offset] = local_max;
            minimos[id / offset] = local_min;

            if (offset == T) {
                maxD = maximos[0];
                minA = minimos[0];
            }

        }
        pthread_barrier_wait(&barrera); //esperar a que se terminen de escribir las variables compartidas
        offset *= 2; //reduce a la mitad los procesos activos
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

    maximos = (double*)malloc(T * sizeof(double));
    minimos = (double*)malloc(T * sizeof(double));

    // Inicializacion
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i * N + j] = i * N + j;
            B[j * N + i] = 1.0; //ordenada por columnas
            AB[i * N + j] = 0.0;
            C[j * N + i] = 1.0; //ordenada por columnas
            ABC[i * N + j] = 0.0;
            D[i * N + j] = i * N + j;
            DC[i * N + j] = 0.0;
            DCB[i * N + j] = 0.0;
        }
    }


    //tomar tiempo start
    timetick = dwalltime();
    pthread_barrier_init(&barrera, NULL, T);
    for (int i = 0; i < T; i++) {
        threads_ids[i] = i;
        pthread_create(&hilos[i], NULL, &multBloques, (void*)&threads_ids[i]);
    }
    for (int i = 0; i < T; i++) {
        pthread_join(hilos[i], NULL);
    }
    pthread_barrier_destroy(&barrera);

    double totalTime = dwalltime() - timetick;
    printf("Tiempo: %fs\n", totalTime);

    // //estoy multiplicando bien?
    // printf("Matriz A\n");
    // printMatriz(A, N);
    // printf("Matriz AB\n");
    // printMatriz(AB, N);
    // printf("Matriz DC\n");
    // printMatriz(DC, N);
    // printf("Matriz ABC\n");
    // printMatriz(ABC, N);
    // printf("Matriz DCB\n");
    // printMatriz(DCB, N);
    printf("MaximoD: %f, MinimoA:%f\n", maxD, minA);
    return 0;
}