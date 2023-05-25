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

void multBloques(double* A, double* B, double* AB, double* C, double* ABC, double* D, double* DC, double* DCB, int N, int BS, double* maxD, double* minA) {
    int id = omp_get_thread_num();
    int T = omp_get_num_threads();
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


    printf("id: %d el maximo local es: %f y el minimo local es: %f\n", id, local_max, local_min);

    // Reducción para obtener el mínimo y maximo 
    #pragma omp critical
    {
        if (local_min < *minA) {
            *minA = local_min;
        }
        if (local_max > *maxD) {
            *maxD = local_max;
        }
    }
}

int main() {
    double* A, * B, * C, * D, * P, * R, * AB, * ABC, * DC, * DCB;
    double maxD = -1, minA = 999;
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
            ABC[i * N + j] = 0.0;
            D[i * N + j] = 1.0;
            DC[i * N + j] = 0.0;
            DCB[i * N + j] = 0.0;
            P[i * N + j] = 0.0;
            R[i * N + j] = 0.0;
        }
    }

    //tomar tiempo start
    double timetick = dwalltime();

	#pragma omp parallel num_threads(4) reduction(max:maxD) reduction(min:minA)
	{
		multBloques(A,B,AB,C,ABC,D,DC,DCB,N,BS,&maxD,&minA);
        printf("el maximo es: %f y el minimo es: %f\n", maxD, minA);

        //#pragma omp barrier

    }

   	double totalTime = dwalltime() - timetick;
    printf("Tiempo en bloques de %d x %d: %f\n", BS, BS, totalTime);


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