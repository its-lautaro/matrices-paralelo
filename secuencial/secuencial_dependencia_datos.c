/*************************************************
* Operaciones en bloque sobre matrices de NxN
* Secuencial, aprovechando la dependencia de datos
*
* * Torres, Gabriel
* * La Vecchia, Lautaro
*************************************************/
#include <stdio.h>
#include<stdlib.h>  
#include <sys/time.h>

double dwalltime() {
    double sec;
    struct timeval tv;

    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

// Valida los resultados para matrices unitarias (donde todos los elementos contienen el mismo resultado)
void validar(double* matriz, double res, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (matriz[i * N + j] != res) {
                printf("Error en %d, %d, valor: %f\n", i, j, matriz[i * N + j]);
                return;
            }
        }
    }
    printf("OK!\n");
}

// Calcula el producto ABC y DCB, obtiene el minimo de A y el maximo de D
void multBloques(double* A, double* B, double* AB, double* C, double* ABC, double* D, double* DC, double* DCB, int N, int bs, double* max, double* min) {
    int I, J, K, i, j, k;
    double* Ablk, * Bblk, * Cblk, * Dblk, * ABblk, * ABCblk, * DCblk, * DCBblk;

    for (I = 0; I < N; I += bs) {
        for (J = 0; J < N; J += bs) {
            ABblk = &AB[I * N + J];
            DCblk = &DC[I * N + J];
            for (K = 0; K < N; K += bs) {
                Ablk = &A[I * N + K];
                Bblk = &B[J * N + K];
                Dblk = &D[I * N + K];
                Cblk = &C[J * N + K];
                for (i = 0; i < bs; i++) {
                    for (j = 0; j < bs; j++) {
                        for (k = 0; k < bs; k++) {
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
        // Una vez que calcule un bloque de AB (y DC), calculo un bloque de ABC (y DCB)
        for (J = 0; J < N; J += bs) {
            ABCblk = &ABC[I * N + J];
            DCBblk = &DCB[I * N + J];
            for (K = 0; K < N; K += bs) {
                ABblk = &AB[I * N + K];
                Cblk = &C[J * N + K];
                DCblk = &DC[I * N + K];
                Bblk = &B[J * N + K];
                for (i = 0; i < bs; i++) {
                    for (j = 0; j < bs; j++) {
                        for (k = 0; k < bs; k++) {
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
// Calcula la matriz P = ABC*maxD + DCB*minA, y la suma total de los elementos de P
double sumBloques_promedio(double* ABC, double* DCB, double* P, int N, int bs, double min, double max) {
    double sumTotal = 0;
    double* ABCblk, * DCBblk, * Pblk;

    for (int I = 0; I < N; I += bs) {
        for (int J = 0; J < N; J += bs) {
            Pblk = &P[I * N + J];
            ABCblk = &ABC[I * N + J];
            DCBblk = &DCB[I * N + J];
            for (int i = 0; i < bs; i++) {
                for (int j = 0; j < bs; j++) {
                    Pblk[i * N + j] = max * ABCblk[i * N + j] + min * DCBblk[i * N + j];
                    sumTotal += Pblk[i * N + j];
                }
            }
        }
    }
    return (sumTotal / (N * N));
}
// Calcula el producto escalar R = PromP * R
void prod_escalar(double* P, double* R, double promP, int N, int bs) {
    double* Pblk, * Rblk;

    for (int I = 0; I < N; I += bs) {
        for (int J = 0; J < N; J += bs) {
            Rblk = &R[I * N + J];
            Pblk = &P[I * N + J];
            for (int i = 0; i < bs; i++) {
                for (int j = 0; j < bs; j++) {
                    Rblk[i * N + j] = promP * Pblk[i * N + j];
                }
            }
        }
    }
}

int main(int argc, char* argv[]) {
    double* A, * B, * C, * D; // Matrices
    double* P, * R, * AB, * ABC, * DC, * DCB; // Matrices resultado
    double maxD = -1, minA = 9999;
    int N, BS;
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

    // Inicio de la operación
    double timetick = dwalltime();

    //ABC = AB*C, DCB = DC*B, MinA, MaxD
    multBloques(A, B, AB, C, ABC, D, DC, DCB, N, BS, &maxD, &minA);
    //P = MaxD*ABC + MinA*DCB, PromP
    promP = sumBloques_promedio(ABC, DCB, P, N, BS, minA, maxD);
    //R=PromP*P
    prod_escalar(P, R, promP, N, BS);

    // Calcular tiempo de ejecucion
    double totalTime = dwalltime() - timetick;
    printf("Multiplicación de matriz %dx%d en bloques de %dx%d\n", N, N, BS, BS);
    printf("Tiempo: %.3fs\n\n", totalTime);

    //Validar los resultados
    printf("Validando matriz AB... ");
    validar(AB, N, N);
    printf("Validando matriz DC... ");
    validar(DC, N, N);
    printf("Validando matriz ABC... ");
    validar(ABC, N * N, N);
    printf("Validando matriz DCB... ");
    validar(DCB, N * N, N);
    printf("Validando matriz P... ");
    validar(P, (N * N * maxD) + (N * N * minA), N);
    printf("Validando matriz R... ");
    validar(R, (((N * N) + (N * N)) * promP), N);
    printf("Maximo D:%f, Minimo A:%f, Promedio:%f\n", maxD, minA, promP);

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