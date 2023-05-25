#include <stdio.h>
#include <stdlib.h>  
#include <sys/time.h>

double dwalltime() {
    double sec;
    struct timeval tv;

    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

void validate(double* matriz, double res, int N){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (matriz[i * N  + j] != res){
                printf("Error en %d, %d, valor: %f\n", i, j, matriz[i * N + j]);
            }
        }
    }
}

void printMatriz(double* matriz, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%f ", matriz[i * N + j]);
        }
        printf("\n");
    }
}

void multBloques(double* A, double* B, double* AB, double* C, double* ABC, double* D, double* DC, double* DCB, int N, int BS, double* max, double* min) {
    double* Ablk, * Bblk, * Cblk, * Dblk, * ABblk, * ABCblk, * DCblk, * DCBblk;

    for (int I = 0; I < N; I += BS) {
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
                        if (Ablk[i * N + j] < *min) *min = Ablk[i * N + j];
                        //MAX(D)
                        if (Dblk[i * N + j] > *max) *max = Dblk[i * N + j];
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
}

double sumBloques_promedio(double* ABC, double* DCB, double* P, int N, int BS, double min, double max) {
    double sumTotal = 0;
    double* ABCblk, * DCBblk, * Pblk;

    for (int I = 0; I < N; I += BS) {
        for (int J = 0; J < N; J += BS) {
            Pblk = &P[I * N + J];
            ABCblk = &ABC[I * N + J];
            DCBblk = &DCB[I * N + J];
            for (int i = 0; i < BS; i++) {
                for (int j = 0; j < BS; j++) {
                    Pblk[i * N + j] = max * ABCblk[i * N + j] + min * DCBblk[i * N + j];
                    sumTotal += Pblk[i * N + j];
                }
            }
        }
    }

    //printf("sumaTotal: %f",sumTotal);
    return (sumTotal / (N * N));
}

void prod_escalar(double* P, double* R, double promP, int N, int BS) {
    double* Pblk, * Rblk;

    for (int I = 0; I < N; I += BS) {
        for (int J = 0; J < N; J += BS) {
            Rblk = &R[I * N + J];
            Pblk = &P[I * N + J];
            for (int i = 0; i < BS; i++) {
                for (int j = 0; j < BS; j++) {
                    Rblk[i * N + j] = promP * Pblk[i * N + j];
                }
            }
        }
    }
}

int main(int argc, char* argv[]) {
    double* A, * B, * C, * D, * P, * R, * AB, * ABC, * DC, * DCB;
    double maxD = 1, minA = 1;
    int N,BS;
    int i, j;
    double promP;

    // Chequeo de parámetros
    if ((argc != 3) || ((N = atoi(argv[1])) <= 0) || ((BS = atoi(argv[2])) <= 0) || ((N % BS) != 0)) {
        printf("Error en los parámetros. Usar: ./%s N BS (N debe ser multiplo de BS)\n", argv[0]);
        exit(1);
    }

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
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
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

    //ABC = AB*C, DCB = DC*B, MinA, MaxD
    //multBloques(A, B, AB, C, ABC, D, DC, DCB, N, BS, &maxD, &minA);
    //P = MaxD*ABC + MinA*DCB, PromP
    promP = sumBloques_promedio(ABC, DCB, P, N, BS, minA, maxD);

    prod_escalar(P,R,promP,N,BS);

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