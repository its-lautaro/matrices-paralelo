#include <stdio.h>
#include <stdlib.h>  
#include <sys/time.h>
#include <omp.h>

double dwalltime() {
    double sec;
    struct timeval tv;

    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

void printMatriz(double* matriz, int N) {
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%f ", matriz[i * N + j]);
        }
        printf("\n");
    }
}


void sum_promedio(double* ABC, double* DCB, double* P, int N, double minA, double maxD, double* promP, int BS) {
    int id = omp_get_thread_num();
    int T = omp_get_num_threads();
    double* Pblk, * ABCblk, * DCBblk;
    double sum_local = 0;

    //for (int I = (id*(N/T)); I < ((id+1)*(N/T)); I ++) { otra opcion 
    for (int I = (id * BS); I < N; I += (T * BS)) {
        for (int J = 0; J < N; J += BS) {
            Pblk = &P[I * N + J];
            ABCblk = &ABC[I * N + J];
            DCBblk = &DCB[I * N + J];
            for (int i = 0; i < BS; i++) {
                for (int j = 0; j < BS; j++) {
                    Pblk[i * N + j] = maxD * ABCblk[i * N + j] + minA * DCBblk[i * N + j];
                    sum_local += Pblk[i * N + j];
                }
            }
        }
    }

    // Reducción para obtener el promedio
    #pragma omp critical
    {
        *promP += sum_local;
    }
}

void producto_escalar(double* P, double* R, int BS, int N, double* promP) {
    int id = omp_get_thread_num();
    int T = omp_get_num_threads();
    double* Pblk, * Rblk;

    for (int I = (id * BS); I < N; I += (T * BS)) {
        for (int J = 0; J < N; J += BS) {
            Rblk = &R[I * N + J];
            Pblk = &P[I * N + J];
            for (int i = 0; i < BS; i++) {
                for (int j = 0; j < BS; j++) {
                    Rblk[i * N + j] = *promP * Pblk[i * N + j];
                }
            }
        }
    }
}

int main() {
    double* A, * B, * C, * D, * P, * R, * AB, * ABC, * DC, * DCB;
    double maxD = 1, minA = 1;
    int BS = 4;
    int N = 4;
    double promP = 0;

    // Alocar  
    A = (double*)malloc(N * N * sizeof(double));
    B = (double*)malloc(N * N * sizeof(double));
    AB = (double*)malloc(N * N * sizeof(double));
    C = (double*)malloc(N * N * sizeof(double));
    ABC = (double*)malloc(N * N * sizeof(double));
    D = (double*)malloc(N * N * sizeof(double));
    DC = (double*)malloc(N * N * sizeof(double));
    DCB = (double*)malloc(N * N * sizeof(double));
    P = (double*)malloc(N * N * sizeof(double));
    R = (double*)malloc(N * N * sizeof(double));

    // Inicializacion
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i * N + j] = 1.0;
            B[j * N + i] = 1.0; //ordenada por columnas
            AB[i * N + j] = 0.0;
            C[j * N + i] = 1.0; //ordenada por columnas
            ABC[i * N + j] = i*N+j;
            D[i * N + j] = 1.0;
            DC[i * N + j] = 0.0;
            DCB[i * N + j] = i*N+j;
            P[i * N + j] = 0.0;
            R[i * N + j] = 0.0;
        }
    }

    //tomar tiempo start
    double timetick = dwalltime();

    #pragma omp parallel num_threads(4) reduction(+:promP) 
    {

        //P = MaxD*ABC + MinA*DCB, PromP
        sum_promedio(ABC,DCB,P,N,minA,maxD,&promP,BS);

        #pragma omp barrier

        #pragma omp single 
        {
            promP = promP/(N*N); 
        }

        producto_escalar(P,R,BS,N,&promP);
    }

    double totalTime = dwalltime() - timetick;
    printf("Tiempo en bloques de %d x %d: %f\n", BS, BS, totalTime);

    printMatriz(R,N);


    //liberamos memoria
    free(A);
    free(B);
    free(C);
    free(D);
    free(AB);
    free(DC);
    free(ABC);
    free(DCB);
    free(P);
    free(R);

    return 0;
}