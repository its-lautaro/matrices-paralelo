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

// void printMatriz(double* matriz, int N) {
//     int i, j;

//     for (i = 0; i < N; i++) {
//         for (j = 0; j < N; j++) {
//             printf("%f ", matriz[i * N + j]);
//         }
//         printf("\n");
//     }
// }

void validarAB(double* matriz, int N) {
    for (int i = 0; i < N; i++) {
        double res = 0; //valor utilizado para validar
        for (int j = 1; j < N + 1; j++) {
            res += i * N + j;
        }
        for (int j = 0; j < N; j++) {
            if (matriz[i * N + j] != res) {
                printf("ERROR %f no coincide con %f\n", matriz[i * N + j], res);
                return;
            }
        }
    }
    printf("OK\n");
}

void validarABC(double* matriz, int N) {
    double res;
    for (int i = 0; i < N; i++) {
        res = 0; //valor utilizado para validar
        for (int j = 0; j < N; j++) {
            res += (i * N + j) + 1; //elemento de AB fila i
        }
        res *= N; // elemento de ABC fila i
        for (int j = 0; j < N; j++) {
            if (matriz[i * N + j] != res) {
                printf("ERROR %f no coincide con %f\n", matriz[i * N + j], res);
                return;
            }
        }
    }
    printf("OK\n");
}

void validarP(double* matriz, int N) {
    double res;
    for (int i = 0; i < N; i++) {
        res = 0; //valor utilizado para validar
        for (int j = 0; j < N; j++) {
            res += (i * N + j) + 1; //elemento de AB fila i
        }
        res *= N; // elemento de ABC fila i
        res = res * 1 + res * (N * N); //elemento de P fila i
        for (int j = 0; j < N; j++) {
            if (matriz[i * N + j] != res) {
                printf("ERROR %f no coincide con %f\n", matriz[i * N + j], res);
                return;
            }
        }
    }
    printf("OK\n");
}

void validarR(double* matriz, int N, double promP) {
    double res;
    for (int i = 0; i < N; i++) {
        res = 0; //valor utilizado para validar
        for (int j = 1; j < N + 1; j++) {
            res += i * N + j;
        }
        res *= N; // elemento de ABC fila i
        res = res * 1 + res * (N * N); //elemento de P fila i
        res = res * promP; //elemento de R fila i
        for (int j = 0; j < N; j++) {
            if (matriz[i * N + j] != res) {
                printf("ERROR %f no coincide con %f\n", matriz[i * N + j], res);
                return;
            }
        }
    }
    printf("OK\n");
}

void multBloques(double* A, double* B, double* AB, double* C, double* ABC, double* D, double* DC, double* DCB, int N, int BS, int id, int T, double* minimos, double* maximos) {
    double* Ablk, * Bblk, * Cblk, * Dblk, * ABblk, * ABCblk, * DCblk, * DCBblk;
    double min_local = 9999, max_local = -1;

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
                        if (Ablk[i * N + j] < min_local) {
                            min_local = Ablk[i * N + j];
                            //printf("proceso %d local_min %f\n", id, local_min);
                        }

                        //MAX(D)
                        if (Dblk[i * N + j] > max_local) {
                            max_local = Dblk[i * N + j];
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
    minimos[id] = min_local;
    maximos[id] = max_local;
}

void sum_promedio(double* ABC, double* DCB, double* P, int N, double minA, double maxD, double* sumas, int BS, int id, int T) {
    double* Pblk, * ABCblk, * DCBblk;
    double sum = 0;

    //for (int I = (id*(N/T)); I < ((id+1)*(N/T)); I ++) { otra opcion 
    for (int I = (id * BS); I < N; I += (T * BS)) {
        for (int J = 0; J < N; J += BS) {
            Pblk = &P[I * N + J];
            ABCblk = &ABC[I * N + J];
            DCBblk = &DCB[I * N + J];
            for (int i = 0; i < BS; i++) {
                for (int j = 0; j < BS; j++) {
                    Pblk[i * N + j] = maxD * ABCblk[i * N + j] + minA * DCBblk[i * N + j];
                    sum += Pblk[i * N + j];
                }
            }
        }
    }

    sumas[id] = sum;
}

void producto_escalar(double* P, double* R, int BS, int id, int T, int N, double promP) {
    double* Pblk, * Rblk;

    for (int I = (id * BS); I < N; I += (T * BS)) {
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

int main() {
    int BS = 64;
    int N = 1024;
    int T = 4;
    double* A, * B, * C, * D;
    double* P, * R, * AB, * ABC, * DC, * DCB;
    double* sumas, * minimos, * maximos;
    double maxD = -1.0, minA = 9999.0;
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
    sumas = (double*)malloc(T * sizeof(double));
    maximos = (double*)malloc(T * sizeof(double));
    minimos = (double*)malloc(T * sizeof(double));

    // Inicializacion
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i * N + j] = (i * N + j) + 1;
            B[j * N + i] = 1.0; //ordenada por columnas
            AB[i * N + j] = 0.0;
            C[j * N + i] = 1.0; //ordenada por columnas
            ABC[i * N + j] = 0.0;
            D[i * N + j] = (i * N + j) + 1;
            DC[i * N + j] = 0.0;
            DCB[i * N + j] = 0.0;
            P[i * N + j] = 0.0;
            R[i * N + j] = 0.0;
        }
    }

    //tomar tiempo start
    double timetick = dwalltime();

    omp_set_num_threads(T);
#pragma omp parallel    
    {
        int id = omp_get_thread_num();
        multBloques(A, B, AB, C, ABC, D, DC, DCB, N, BS, id, T, minimos, maximos);

        //Reducción para obtener el mínimo y maximo
#pragma omp critical 
        {
            if (minimos[id] < minA) {
                minA = minimos[id];
            }
            if (maximos[id] > maxD) {
                maxD = maximos[id];
            }
        }

#pragma omp barrier

        //P = MaxD*ABC + MinA*DCB, PromP
        sum_promedio(ABC, DCB, P, N, minA, maxD, sumas, BS, id, T);

        //Reducción para obtener el promedio
#pragma omp critical
        {
            promP += sumas[id];
        }

#pragma omp barrier

#pragma omp single 
        {
            promP = promP / (N * N);
        }

        producto_escalar(P, R, BS, id, T, N, promP);
    }

    double totalTime = dwalltime() - timetick;
    printf("Tiempo en bloques de %d x %d: %f\n", BS, BS, totalTime);
    printf("Maximo %f Minimo %f PromP %f\n", maxD, minA, promP);
    printf("Validando producto AB.. ");
    validarAB(AB, N);
    printf("Validando producto DC.. ");
    validarAB(DC, N);
    printf("Validando producto ABC.. ");
    validarABC(ABC, N);
    printf("Validando producto DCB.. ");
    validarABC(DCB, N);
    //maximo (N*N) minimo 1
    printf("Validando matriz P.. ");
    validarP(P, N);
    printf("Validando matriz R.. ");
    validarR(R, N, promP);

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
    free(sumas);
    free(maximos);
    free(minimos);

    return 0;
}