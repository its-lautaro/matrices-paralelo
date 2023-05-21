#include <stdio.h>

int maxMatriz(int matriz[][], int n) {
    int maxVal = matriz[0][0];
    int i, j;
    
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (matriz[i][j] > maxVal) {
                maxVal = matriz[i][j];
            }
        }
    }
    
    return maxVal;
}

int minMatriz(int matriz[][], int n) {
    int minVal = matriz[0][0];
    int i, j;
    
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (matriz[i][j] < minVal) {
                minVal = matriz[i][j];
            }
        }
    }
    
    return minVal;
}

void multBloques(int A[][], int B[][], int C[][], int n, int bs) {
    int i, j, k, ibloque, jbloque, kbloque;
    
    for (ibloque = 0; ibloque < n; ibloque += bs) {
        for (jbloque = 0; jbloque < n; jbloque += bs) {
            for (kbloque = 0; kbloque < n; kbloque += bs) {
                for (i = ibloque; i < ibloque + bs; i++) {
                    for (j = jbloque; j < jbloque + bs; j++) {
                        for (k = kbloque; k < kbloque + bs; k++) {
                            C[i][j] += A[i][k] * B[k][j];
                        }
                    }
                }
            }
        }
    }
}

void sumBloques(int A[][], int B[][], int C[][], int n, int bs) {
    int i, j, ibloque, jbloque;
    
    for (ibloque = 0; ibloque < n; ibloque += bs) {
        for (jbloque = 0; jbloque < n; jbloque += bs) {
            for (i = ibloque; i < ibloque + bs; i++) {
                for (j = jbloque; j < jbloque + bs; j++) {
                    C[i][j] = A[i][j] + B[i][j];
                }
            }
        }
    }
}

void multEscalar(int matriz[][], int n, int escalar, int bs) {
    int i, j, bi, bj;

    for (bi = 0; bi < n; bi += bs) {
        for (bj = 0; bj < n; bj += bs) {
            for (i = bi; i < bi + bs && i < n; i++) {
                for (j = bj; j < bj + bs && j < n; j++) {
                    matriz[i][j] *= escalar;
                }
            }
        }
    }
}

void multProm(int matriz[][], int result[][], int n, float prom, int bs) {
    int i, j, bi, bj;

    for (bi = 0; bi < n; bi += bs) {
        for (bj = 0; bj < n; bj += bs) {
            for (i = bi; i < bi + bs && i < n; i++) {
                for (j = bj; j < bj + bs && j < n; j++) {
                    result[i][j] = matriz[i][j]*prom;
                }
            }
        }
    }
}

float prom(int m[][], int n) {
    int i, j;
    int sum = 0;
    int total = n * n;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            sum += m[i][j];
        }
    }

    return (float)sum/total;
}

void printMatriz(int matriz[][], int n) {
    int i, j;
    
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%d ", matriz[i][j]);
        }
        printf("\n");
    }
}

int main() {
    double *A, *B, *C, *D, *P, *R, *AB, *ABC, *DC, *DCB;
    int maxD, minA;
    int n, bs;
    int i,j,k;
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
    ABC = (double*)malloc(n * n * sizeof(double));
    C = (double*)malloc(n * n * sizeof(double));
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
    
    int result[][];
    
    multBloques(A, B, AB, n, bs); //A*B
    multBloques(AB, C, ABC, n, bs); //AB*C

    multBloques(D, C, DC, n, bs); //D*C
    multBloques(DC, B, DCB, n, bs); //DC*B

    maxD = maxMatriz(D, n);
    minA = minMatriz(A, n);
    
    //multiplicamos las matrices por los valores maximos y minimos respectivos
    multEscalar(ABC, n, maxD, bs);
    multEscalar(DCB, n, minD, bs);

    //realizamos la suma y lo guardamos en P
    sumBloques(ABC, DCB, P, n, bs);

    //printf("El resultado de la operación P = MaxD.(ABC) + MinA.(DCB) es:\n");
    //printMatriz(P, n);

    promP = prom(P, n);

    //multiplicamos a P por su prom
    multProm(P, R, promP, bs);

    //printf("El resultado de la operación R = Prom(P)*P es:\n");
    //printMatriz(R, n);


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
