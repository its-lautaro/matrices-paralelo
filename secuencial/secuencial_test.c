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

// Calcula el producto matricial de dos matrices de tamaño nxn, en bloques de tamaño bsxbs
void multBloques(double* A, double* B, double* C, int n, int bs) {
    double* ablk, * bblk, * cblk;
    for (int I = 0; I < n; I += bs) {
        for (int J = 0; J < n; J += bs) {
            cblk = &C[I * n + J];
            for (int K = 0; K < n; K += bs) {
                ablk = &A[I * n + K];
                bblk = &B[J * n + K];
                for (int i = 0; i < bs; i++) {
                    for (int j = 0; j < bs; j++) {
                        for (int k = 0; k < bs; k++) {
                            cblk[i * n + j] += ablk[i * n + k] * bblk[j * n + k];
                        }
                    }
                }
            }
        }
    }
}

// Suma dos matrices cuadradas, en bloques
void sumBloques(double* A, double* B, double* C, int n, int bs) {
    double* ablk, * bblk, * cblk;
    for (int I = 0; I < n; I += bs) {
        for (int J = 0; J < n; J += bs) {
            cblk = &C[I * n + J];
            ablk = &A[I * n + J];
            bblk = &B[I * n + J];
            for (int i = 0; i < bs; i++) {
                for (int j = 0; j < bs; j++) {
                    cblk[i * n + j] = ablk[i * n + j] + bblk[i * n + j];
                }
            }
        }
    }
}

// Realiza el producto escalar a la matriz, en bloques
void multEscalar(double* matriz, int n, int escalar, int bs) {
    int i, j, bi, bj;
    //multiplicacion en bloques
    for (bi = 0; bi < n; bi += bs) {
        for (bj = 0; bj < n; bj += bs) {
            for (i = bi;(i < bi + bs) && (i < n); i++) {
                for (j = bj; (j < bj + bs) && (j < n); j++) {
                    matriz[i * n + j] *= escalar;
                }
            }
        }
    }
}


int main(int argc, char* argv[]) {
    double* A, * B, * C, * D; // Matrices
    double* P, * AB, * ABC, * DC, * DCB; // Matrices resultado
    double maxD = -1, minA = 9999, sum = 0;
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
        }
    }

    // Inicio de la operación
    double timetick = dwalltime();

    // Realizamos la multiplicacion de matrices
    multBloques(A, B, AB, N, BS);
    multBloques(AB, C, ABC, N, BS);
    multBloques(D, C, DC, N, BS);
    multBloques(DC, B, DCB, N, BS);
    // Calculamos el valor maximo de D
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (D[i * N + j] > maxD) {
                maxD = D[i * N + j];
            }
        }
    }

    // Calculamos el valor minimo de A
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (A[i * N + j] < minA) {
                minA = A[i * N + j];
            }
        }
    }

    // Multiplicamos las matrices por los valores maximos y minimos respectivos
    multEscalar(ABC, N, maxD, BS); //ABC*maxD = (n*n*max)
    multEscalar(DCB, N, minA, BS); //DCB*miniD = (n*n*min)

    // Realizamos la suma y lo guardamos en P
    sumBloques(ABC, DCB, P, N, BS); //minA*ABC + maxD*DCB = (n*n*2)

    // Calculamos el promedio de P
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            sum += P[i * N + j];
        }
    }

    //multiplicamos a P por su promedio
    promP = (sum / (N * N));
    multEscalar(P, N, promP, BS);

    //Calcular tiempo de ejecucion
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
    validar(P, (((N * N) + (N * N)) * promP), N);
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

    return 0;
}
