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

double maxMatriz(double* matriz, int n) {
    double maxVal = -1;
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (matriz[i * n + j] > maxVal) {
                maxVal = matriz[i * n + j];
            }
        }
    }

    return maxVal;
}

double minMatriz(double* matriz, int n) {
    double minVal = 9999;
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (matriz[i * n + j] < minVal) {
                minVal = matriz[i * n + j];
            }
        }
    }

    return minVal;
}

void multBloques(double* A, double* B, double* C, int n, int bs) {
    int I, J, K, i, j, k;
    double* ablk, * bblk, * cblk;
    for (I = 0; I < n; I += bs) {
        for (J = 0; J < n; J += bs) {
            cblk = &C[I * n + J];
            for (K = 0; K < n; K += bs) {
                ablk = &A[I * n + K];
                bblk = &B[J * n + K];
                for (i = 0; i < bs; i++) {
                    for (j = 0; j < bs; j++) {
                        for (k = 0; k < bs; k++) {
                            cblk[i * n + j] += ablk[i * n + k] * bblk[j * n + k];
                        }
                    }
                }
            }
        }
    }
}

void sumBloques(double* A, double* B, double* C, int n, int bs) {
    //suma en bloques
    int I, J, K, i, j, k;
    double* ablk, * bblk, * cblk;
    for (I = 0; I < n; I += bs) {
        for (J = 0; J < n; J += bs) {
            cblk = &C[I * n + J];
            for (K = 0; K < n; K += bs) {
                ablk = &A[I * n + K];
                bblk = &B[J * n + K];
                for (i = 0; i < bs; i++) {
                    for (j = 0; j < bs; j++) {
                        for (k = 0; k < bs; k++) {
                            cblk[i * n + j] += ablk[i * n + k] + bblk[i * n + k];
                        }
                    }
                }
            }
        }
    }
    //suma
    // for (int i = 0; i < n; i++)
    // {
    //     for (int j = 0; j < n; j++)
    //     {
    //         C[i*n+j] += A[i*n+j] + B[i*n+j];
    //     }

    // }   
}

void multEscalar(double* matriz, int n, int escalar, int bs) {
    int i, j, bi, bj;
    //multiplicacion en bloques
    for (bi = 0; bi < n - bs; bi += bs) {
        for (bj = 0; bj < n - bs; bj += bs) {
            for (i = bi; i < bi + bs; i++) {
                for (j = bj; j < bj + bs; j++) {
                    matriz[i * n + j] *= escalar;
                }
            }
        }
    }

    // multiplicacion
    // for (int i = 0; i < n; i++)
    // {
    //     for (int j = 0; j < n; j++)
    //     {
    //         matriz[i*n+j] *= escalar;
    //     }

    // }

}

void multProm(double* matriz, double* result, int n, int prom, int bs) {
    int i, j, bi, bj;

    for (bi = 0; bi < n - bs; bi += bs) {
        for (bj = 0; bj < n - bs; bj += bs) {
            for (i = bi; i < bi + bs; i++) {
                for (j = bj; j < bj + bs; j++) {
                    result[i * n + j] = matriz[i * n + j] * prom;
                }
            }
        }
    }
}

double prom(double* m, int n) {
    int i, j;
    double sum = 0;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            sum += m[i * n + j];
        }
    }

    return (sum / n * n);
}

void printMatriz(double* matriz, int n) {
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%f ", matriz[i * n + j]);
        }
        printf("\n");
    }
}

int main(int argc, char* argv[]) {
    double* A, * B, * C, * D, * P, * R, * AB, * ABC, * DC, * DCB;
    double maxD, minA;
    int n, bs;
    int i, j, k;
    float promP;

    // Chequeo de parámetros
    if ((argc != 3) || ((n = atoi(argv[1])) <= 0) || ((bs = atoi(argv[2])) <= 0) || ((n % bs) != 0)) {
        printf("Error en los parámetros. Usar: ./%s N BS (N debe ser multiplo de BS)\n", argv[0]);
        exit(1);
    }

    // Alocar  
    A = (double*)malloc(n * n * sizeof(double));
    B = (double*)malloc(n * n * sizeof(double));
    AB = (double*)malloc(n * n * sizeof(double));
    C = (double*)malloc(n * n * sizeof(double));
    ABC = (double*)malloc(n * n * sizeof(double));
    D = (double*)malloc(n * n * sizeof(double));
    DC = (double*)malloc(n * n * sizeof(double));
    DCB = (double*)malloc(n * n * sizeof(double));
    P = (double*)malloc(n * n * sizeof(double));
    R = (double*)malloc(n * n * sizeof(double));

    // Inicializacion
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            A[i * n + j] = 1.0;
            B[j * n + i] = 1.0; //ordenada por columnas
            AB[i * n + j] = 0.0;
            C[j * n + i] = 1.0; //ordenada por columnas
            ABC[i * n + j] = 0.0;
            D[i * n + j] = 1.0;
            DC[i * n + j] = 0.0;
            DCB[i * n + j] = 0.0;
            P[i * n + j] = 0.0;
            R[i * n + j] = 0.0;
        }
    }
    //tomar tiempo start
    //double timetick = dwalltime();

    multBloques(A, B, AB, n, bs); //A*B (todos los elementos n)
    multBloques(AB, C, ABC, n, bs); //AB*C 

    multBloques(D, C, DC, n, bs); //D*C
    multBloques(DC, B, DCB, n, bs); //DC*B

    maxD = maxMatriz(D, n);
    minA = minMatriz(A, n);
    printf("El maximo es %d el minimo es %d\n", maxD, minA);

    //multiplicamos las matrices por los valores maximos y minimos respectivos
    multEscalar(ABC, n, maxD, bs);
    printf("El resultado de ABC es:\n");
    printMatriz(ABC, n);
    
    multEscalar(DCB, n, minA, bs);
    printf("El resultado de DCB es:\n");
    printMatriz(DCB, n);
    //realizamos la suma y lo guardamos en P
    sumBloques(ABC, DCB, P, n, bs);

    //double totalTime = dwalltime() - timetick;
    //printf("Tiempo en bloques de %d x %d: %f\n", bs, bs, totalTime);

    printf("El resultado de la operación P = MaxD.(ABC) + MinA.(DCB) es:\n");
    printMatriz(P, n);

    promP = prom(P, n);

    //consultar
    //multiplicamos a P por su prom
    multProm(P, R, n, promP, bs);
    //multEscalar(P, n, promP, bs);

    printf("El resultado de la operación R = Prom(P)*P es:\n");
    printMatriz(R, n);

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
