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

void validate(double* matriz, double res, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (matriz[i * N + j] != res) {
                printf("Error en %d, %d, valor: %f\n", i, j, matriz[i * N + j]);
            }
        }
    }
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
            ablk = &A[I * n + J];
            bblk = &B[I * n + J];
            for (i = 0; i < bs; i++) {
                for (j = 0; j < bs; j++) {
                    cblk[i * n + j] = ablk[i * n + j] + bblk[i * n + j];
                }
            }
        }
    }
}

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
    double* A, * B, * C, * D, * P, * AB, * ABC, * DC, * DCB;
    double maxD = -1, minA = 9999, sum = 0;
    int n, bs;
    int i, j, k;
    double promP;

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
    //R = (double*)malloc(n * n * sizeof(double));

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
            //R[i * n + j] = 0.0;
        }
    }

    //tomar tiempo start
    double timetick = dwalltime();

    //realizamos la multiplicacion de matrices
    multBloques(A, B, AB, n, bs); //A*B (n)
    multBloques(AB, C, ABC, n, bs); //AB*C (n*n)
    multBloques(D, C, DC, n, bs); //D*C (n)
    multBloques(DC, B, DCB, n, bs); //DC*B (n*n)

    //calculamos el valor maximo de D
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (D[i * n + j] > maxD) {
                maxD = D[i * n + j];
            }
        }
    }

    //calculamos el valor minimo de A
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (A[i * n + j] < minA) {
                minA = A[i * n + j];
            }
        }
    }

    //multiplicamos las matrices por los valores maximos y minimos respectivos
    multEscalar(ABC, n, maxD, bs); //ABC*maxD = (n*n*max)
    multEscalar(DCB, n, minA, bs); //DCB*miniD = (n*n*min)

    //realizamos la suma y lo guardamos en P
    sumBloques(ABC, DCB, P, n, bs); //minA*ABC + maxD*DCB = (n*n*2)

    //calculamos el promedio de P
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sum += P[i * n + j];
        }
    }

    //multiplicamos a P por su promedio
    promP = (sum / (n * n));
    multEscalar(P, n, promP, bs);

    double totalTime = dwalltime() - timetick;
    printf("Tiempo en bloques de %d x %d: %f\n", bs, bs, totalTime);

    //Validaciones 
    //AB = A*B , todos sus valores son N
    validate(AB, n, n);
    //DC = D*C , todos sus valores son N
    validate(DC, n, n);
    //ABC = AB*C , todos sus valores son N*N*maxD
    validate(ABC, n * n * maxD, n);
    //DCB= DC*B , todos sus valores son N*N*minA
    validate(DCB, n * n * minA, n);
    //R = promP*P, todos sus valores son promP*((N*N*maxD)+(N*N*minA))
    validate(P, ((n * n * maxD) + (n * n * minA)) * (sum / (n * n)), n);
    
    //validate(R, (promP * ((n * n * maxD) + (n * n * minA))), n);

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
    //free(R);

    return 0;
}
