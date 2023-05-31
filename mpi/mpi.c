/*************************************************
* Operaciones en bloque sobre matrices de NxN
* Comunicación colectiva en MPI
*
* * Torres, Gabriel
* * La Vecchia, Lautaro
*************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>

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
void producto_matricial(double* A, double* B, double* AB, double* D, double* C, double* DC, double* ABC, double* DCB, double* max, double* min, int N, int BS, int filas) {
    double local_min = N * N;
    double local_max = -1;

    for (int i = 0; i < filas; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                AB[i * N + j] += A[i * N + k] * B[j * N + k];
                DC[i * N + j] += D[i * N + k] * C[j * N + k];
            }
            //MIN(A)
            if (A[i * N + j] < local_min) {
                local_min = A[i * N + j];
            }

            //MAX(D)
            if (D[i * N + j] > local_max) {
                local_max = D[i * N + j];
            }
        }
        // Una vez que calcule un bloque de AB (y DC), calculo un bloque de ABC (y DCB)
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                ABC[i * N + j] += AB[i * N + k] * C[j * N + k];
                DCB[i * N + j] += DC[i * N + k] * B[j * N + k];
            }
        }
    }
    MPI_Reduce(&local_max, max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_min, min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
}

// Calcula la matriz P = ABC*maxD + DCB*minA, y la suma total de los elementos de P
double sum_promedio(double* ABC, double* DCB, double* P, int N, double minA, double maxD, double* sumTotal, int BS, int filas) {
    double sum = 0;

    for (int I = 0; I < filas; I += BS) {
        for (int J = 0; J < N; J += BS) {
            for (int i = I; i < I + BS; i++) {
                for (int j = J; j < J + BS; j++) {
                    P[i * N + j] = maxD * ABC[i * N + j] + minA * DCB[i * N + j];
                    sum += P[i * N + j];
                }
            }
        }
    }
    // Realizo la reduccion de las sumas
    MPI_Reduce(&sum, sumTotal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

// Calcula el producto escalar R = PromP * R
void producto_escalar(double* P, double* R, int N, double promP, int BS, int filas) {
    for (int I = 0; I < filas; I += BS) {
        for (int J = 0; J < N; J += BS) {
            for (int i = I; i < I + BS; i++) {
                for (int j = J; j < J + BS; j++) {
                    R[i * N + j] = promP * P[i * N + j];
                }
            }
        }
    }
}

int main(int argc, char** argv) {
    int id;
    int nro_procesos;
    int N; // Dimension
    int BS; // Tamaño de bloque
    double timetick;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &nro_procesos);

    // Chequeo de parametros
    if ((argc != 3) || ((N = atoi(argv[1])) <= 0) || ((BS = atoi(argv[2])) <= 0) || ((N % BS) != 0)) {
        printf("Error en los parámetros. Usar: %s N BS (N debe ser multiplo de BS)\n", argv[0]);
        return -1;
    }

    int filas_por_proceso = N / nro_procesos;
    int elementos_por_proceso = filas_por_proceso * N;
    double* A, * B, * C, * D; // Matrices
    double* R, * P, * AB, * DC, * DCB, * ABC; // Matrices resultado
    double maxD, minA;
    double sumTotal;
    double promP;
    // Bloques
    double* Ablk = (double*)malloc(elementos_por_proceso * sizeof(double));
    double* Cblk = (double*)malloc(elementos_por_proceso * sizeof(double));
    double* Dblk = (double*)malloc(elementos_por_proceso * sizeof(double));
    double* Pblk = (double*)malloc(elementos_por_proceso * sizeof(double));
    double* Rblk = (double*)malloc(elementos_por_proceso * sizeof(double));
    // Bloques resultado (inicializados en 0)
    double* ABblk = (double*)calloc(elementos_por_proceso, sizeof(double));
    double* DCblk = (double*)calloc(elementos_por_proceso, sizeof(double));
    double* DCBblk = (double*)calloc(elementos_por_proceso, sizeof(double));
    double* ABCblk = (double*)calloc(elementos_por_proceso, sizeof(double));
    // Matrices enteras (para el producto)
    B = (double*)malloc(N * N * sizeof(double)); // Todos los procesos necesitan la matriz B
    C = (double*)malloc(N * N * sizeof(double)); // Todos los procesos necesitan la matriz C

    if (id == 0) {
        // Alocar
        A = (double*)malloc(N * N * sizeof(double));
        D = (double*)malloc(N * N * sizeof(double));

        // Inicializar
        for (int i = 0; i < N * N; i++) {
            A[i] = 1.0;
            D[i] = 1.0;
            B[i] = 1.0;
            C[i] = 1.0;
        }
    }

    // El proceso root evalua el tiempo de ejecución
    if (id == 0) {
        timetick = dwalltime();
    }
    MPI_Barrier(MPI_COMM_WORLD); // Necesaria para el calculo correcto del tiempo de ejecución
    // Todos los procesos necesitan la matriz B y C entera para el calculo del producto
    MPI_Bcast(B, N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(C, N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Scatter(A, elementos_por_proceso, MPI_DOUBLE, Ablk, elementos_por_proceso, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(D, elementos_por_proceso, MPI_DOUBLE, Dblk, elementos_por_proceso, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    producto_matricial(Ablk, B, ABblk, Dblk, C, DCblk, ABCblk, DCBblk, &maxD, &minA, N, BS, filas_por_proceso);

    MPI_Bcast(&minA, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&maxD, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    sum_promedio(ABCblk, DCBblk, Pblk, N, minA, maxD, &sumTotal, BS, filas_por_proceso);
    // El proceso root calcula el promedio y lo comparte por broadcast
    if (id == 0) promP = sumTotal / (N * N);
    MPI_Bcast(&promP, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // R=promP*P
    producto_escalar(Pblk, Rblk, N, promP, BS, filas_por_proceso);

    MPI_Barrier(MPI_COMM_WORLD); // Necesaria para el calculo correcto del tiempo de ejecución
    // Calcular tiempo de ejecucion
    if (id == 0) {
        double totalTime = dwalltime() - timetick;
        printf("Multiplicación de matriz %dx%d en bloques de %dx%d\n", N, N, BS, BS);
        printf("Tiempo: %.3fs\n\n", totalTime);

        //necesito alocar espacio para validar las matrices resultado.
        ABC = (double*)malloc(N * N * sizeof(double));
        DCB = (double*)malloc(N * N * sizeof(double));
        AB = (double*)malloc(N * N * sizeof(double));
        DC = (double*)malloc(N * N * sizeof(double));
        P = (double*)malloc(N * N * sizeof(double));
        R = (double*)malloc(N * N * sizeof(double));
    }
    
    //Recopilo todos los elementos de las matrices resultado para validar
    MPI_Gather(Rblk, elementos_por_proceso, MPI_DOUBLE, R, elementos_por_proceso, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(Pblk, elementos_por_proceso, MPI_DOUBLE, P, elementos_por_proceso, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(ABblk, elementos_por_proceso, MPI_DOUBLE, AB, elementos_por_proceso, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(DCblk, elementos_por_proceso, MPI_DOUBLE, DC, elementos_por_proceso, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(ABCblk, elementos_por_proceso, MPI_DOUBLE, ABC, elementos_por_proceso, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(DCBblk, elementos_por_proceso, MPI_DOUBLE, DCB, elementos_por_proceso, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (id == 0) {
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
    }

    MPI_Finalize();
    return 0;
}