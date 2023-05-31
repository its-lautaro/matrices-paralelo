/*************************************************
* Operaciones en bloque sobre matrices de NxN
* OpenMP
*
* * Torres, Gabriel
* * La Vecchia, Lautaro
*************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
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

int main(int argc, char** argv) {
    int N; // Dimension
    int BS; // Tamaño de bloque
    int T; // Número de hilos a utilizar

    // Chequeo de parametros
    if ((argc != 4) || ((N = atoi(argv[1])) <= 0) || ((BS = atoi(argv[2])) <= 0) || ((N % BS) != 0)) {
        printf("Error en los parámetros. Usar: %s N BS T (N debe ser multiplo de BS)\n", argv[0]);
        return -1;
    }
    T=atoi(argv[3]);

    double* A, * B, * C, * D; // Matrices
    double* R, * P, * AB, * DC, * DCB, * ABC; // Matrices resultado
    double maxD = -1, minA = 9999;
    double sumaTotal = 0;
    double promP;

    double maximos[T];
    double minimos[T];
    double sumas[T];

    // Alocación de matrices NxN
    A = (double*)malloc(N * N * sizeof(double));
    B = (double*)malloc(N * N * sizeof(double));
    C = (double*)malloc(N * N * sizeof(double));
    AB = (double*)malloc(N * N * sizeof(double));
    DC = (double*)malloc(N * N * sizeof(double));
    ABC = (double*)malloc(N * N * sizeof(double));
    D = (double*)malloc(N * N * sizeof(double));
    DCB = (double*)malloc(N * N * sizeof(double));
    P = (double*)malloc(N * N * sizeof(double));
    R = (double*)malloc(N * N * sizeof(double));

    // Inicializar las matrices A, B y AB
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i * N + j] = 1.0;
            B[j * N + i] = 1.0;
            AB[i * N + j] = 0;
            C[j * N + i] = 1.0;
            ABC[i * N + j] = 0;
            D[i * N + j] = 1.0;
            DC[i * N + j] = 0;
            DCB[i * N + j] = 0;
        }
    }

    double timetick = dwalltime();
    // Inicio de operación paralela
    #pragma omp parallel num_threads(T)
    {
        // Realizar la multiplicación de matrices en bloques (ABC y DCB) y construir arreglo de maximos y minimos
        #pragma omp for reduction(max: maxD) reduction(min: minA)
        for (int I = 0; I < N; I += BS) {
            for (int J = 0; J < N; J += BS) {
                for (int K = 0; K < N; K += BS) {
                    // Calcular la multiplicación de bloques
                    for (int i = I; i < I + BS; i++) {
                        for (int j = J; j < J + BS; j++) {
                            for (int k = K; k < K + BS; k++) {
                                AB[i * N + j] += A[i * N + k] * B[j * N + k];
                                DC[i * N + j] += D[i * N + k] * C[j * N + k];
                            }
                            // MIN(A)
                            if (A[i * N + j] < minA) {
                                minA = A[i * N + j];
                            }
                            // MAX(D)
                            if (D[i * N + j] > maxD) {
                                maxD = D[i * N + j];
                            }
                        }
                    }
                }
            }
            // Una vez que calcule un bloque de AB (y DC), calculo un bloque de ABC (y DCB)
            for (int J = 0; J < N; J += BS) {
                for (int K = 0; K < N; K += BS) {
                    // Calcular la multiplicación de bloques
                    for (int i = I; i < I + BS; i++) {
                        for (int j = J; j < J + BS; j++) {
                            for (int k = K; k < K + BS; k++) {
                                ABC[i * N + j] += AB[i * N + k] * C[j * N + k];
                                DCB[i * N + j] += DC[i * N + k] * B[j * N + k];
                            }
                        }
                    }
                }
            }
        }
        // Calcular la suma P=maxD*ABC + minA*DCB, en bloques. Calcular la suma total de elementos de P
        #pragma omp for reduction(+: sumaTotal)
        for (int I = 0; I < N; I += BS) {
            for (int J = 0; J < N; J += BS) {
                // Calcular la suma en bloques
                for (int i = I; i < I + BS; i++) {
                    for (int j = J; j < J + BS; j++) {
                        P[i * N + j] = maxD * ABC[i * N + j] + minA * DCB[i * N + j];
                        sumaTotal += P[i * N + j];
                    }
                }
            }
        }

        // Para el producto escalar todos los procesos necesitan promP
        promP = sumaTotal / (N * N);

        // Calcula el producto escalar R=PromP*P
        #pragma omp for
        for (int I = 0; I < N; I += BS) {
            for (int J = 0; J < N; J += BS) {
                // Calcular la suma en bloques
                for (int i = I; i < I + BS; i++) {
                    for (int j = J; j < J + BS; j++) {
                        R[i * N + j] = promP * P[i * N + j];
                    }
                }
            }
        }
    }
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

    // Liberar memoria
    free(A);
    free(B);
    free(AB);
    free(D);
    free(C);
    free(DC);
    free(ABC);
    free(DCB);
    free(P);
    free(R);

    return 0;
}
