/*************************************************
* Operaciones en bloque sobre matrices de NxN
* Pthreads
*
* * Torres, Gabriel
* * La Vecchia, Lautaro
*************************************************/
#include <stdio.h>
#include <stdlib.h>  
#include <sys/time.h>
#include <pthread.h>

// Variables compartidas
int N; // Dimension
int BS; // Tamaño de bloque
int T; // Número de hilos a utilizar
double* A, * B, * C, * D; // Matrices
double* P, * ABC, * DCB, * DC, * AB, * R; // Matrices resultado
double* maximos, * minimos, * sumas;
double maxD, minA;
double promP;

pthread_barrier_t barrera;

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
void multBloques(int id) {
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
                        }

                        //MAX(D)
                        if (Dblk[i * N + j] > local_max) {
                            local_max = Dblk[i * N + j];
                        }
                    }
                }
            }
        }
        // Una vez que calcule un bloque de AB (y DC), calculo un bloque de ABC (y DCB)
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
    // Cargo minimo y maximo local
    minimos[id] = local_min;
    maximos[id] = local_max;
    pthread_barrier_wait(&barrera); // Para esperar a que se terminen de escribir las variables compartidas
    // Reducción de los arreglos (max y min)
    int procesos_activos = T / 2;
    while (procesos_activos > 0) {
        if (id < procesos_activos) {
            local_min = (minimos[id * 2] > minimos[(id * 2) + 1]) ? minimos[(id * 2) + 1] : minimos[id * 2];
            local_max = (maximos[id * 2] < maximos[(id * 2) + 1]) ? maximos[(id * 2) + 1] : maximos[id * 2];
        }
        pthread_barrier_wait(&barrera); // Para esperar a que se terminen de leer las variables compartidas
        if (id < procesos_activos) {
            // Reducción
            maximos[id] = local_max;
            minimos[id] = local_min;
        }
        procesos_activos /= 2;
        pthread_barrier_wait(&barrera); // Para esperar a que se terminen de escribir las variables compartidas
    }
    // El thread 0 actua como "root" y comparte los datos a maxD y minA
    if (id == 0) {
        maxD = maximos[0];
        minA = minimos[0];
    }
}
// Calcula la matriz P = ABC*maxD + DCB*minA, y la suma total de los elementos de P
void sum_promedio(int id) {
    double* Pblk, * ABCblk, * DCBblk, * Ablk, * Bblk, * Cblk, * Dblk, * ABblk, * DCblk;
    double sum_local;
    for (int I = (id * BS); I < N; I += (T * BS)) {
        for (int J = 0; J < N; J += BS) {
            Pblk = &P[I * N + J];
            ABCblk = &ABC[I * N + J];
            DCBblk = &DCB[I * N + J];
            for (int i = 0; i < BS; i++) {
                for (int j = 0; j < BS; j++) {
                    Pblk[i * N + j] = maxD * ABCblk[i * N + j] + minA * DCBblk[i * N + j];
                    sum_local += Pblk[i * N + j];
                }
            }
        }
    }
    // Cada proceso guarda su suma local
    sumas[id] = sum_local;
    pthread_barrier_wait(&barrera); // Para esperar a que se terminen de escribir las variables compartidas
    // Reducción del arreglo (suma)
    int procesos_activos = T / 2;
    while (procesos_activos > 0) {
        if (id < procesos_activos) {
            sum_local = 0;
            sum_local = sumas[id * 2] + sumas[id * 2 + 1];
        }
        pthread_barrier_wait(&barrera); // Para esperar a que se terminen de leer las variables compartidas
        if (id < procesos_activos) {
            sumas[id] = sum_local;
        }
        procesos_activos /= 2;
        pthread_barrier_wait(&barrera); // Para esperar a que se terminen de escribir las variables compartidas
    }
    // El thread 0 actua como "root" y comparte el dato promP
    if (id == 0) {
        promP = sumas[0] / (N * N);
    }
}
// Calcula el producto escalar R = PromP * R
void producto_escalar(int id) {
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

void* behavior(void* arg) {
    int id = *(int*)arg;

    multBloques(id);

    pthread_barrier_wait(&barrera); // Para calcular la suma necesito los valores maxD y minA

    sum_promedio(id);

    pthread_barrier_wait(&barrera); // Para calcular el producto escalar necesito el valor de promP

    producto_escalar(id);

    pthread_exit(NULL);
}

int main(int argc, char* argv[]) {

    // Chequeo de parametros
    if ((argc != 4) || ((N = atoi(argv[1])) <= 0) || ((BS = atoi(argv[2])) <= 0) || ((N % BS) != 0)) {
        printf("Error en los parámetros. Usar: %s N BS T (N debe ser multiplo de BS)\n", argv[0]);
        return -1;
    }
    T = atoi(argv[3]);

    pthread_t hilos[T];
    int threads_ids[T];

    // Alocar  
    A = (double*)malloc(N * N * sizeof(double));
    B = (double*)malloc(N * N * sizeof(double));
    AB = (double*)malloc(N * N * sizeof(double));
    C = (double*)malloc(N * N * sizeof(double));
    D = (double*)malloc(N * N * sizeof(double));
    DC = (double*)malloc(N * N * sizeof(double));
    ABC = (double*)malloc(N * N * sizeof(double));
    DCB = (double*)malloc(N * N * sizeof(double));
    P = (double*)malloc(N * N * sizeof(double));
    R = (double*)malloc(N * N * sizeof(double));
    sumas = (double*)malloc(T * sizeof(double));
    maximos = (double*)malloc(T * sizeof(double));
    minimos = (double*)malloc(T * sizeof(double));

    // Inicializacion
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i * N + j] = 1.0;
            B[j * N + i] = 1.0; //ordenada por columnas
            AB[i * N + j] = 0.0;
            C[j * N + i] = 1.0; //ordenada por columnas
            D[i * N + j] = 1.0;
            DC[i * N + j] = 0.0;
            ABC[i * N + j] = 0.0;
            DCB[i * N + j] = 0.0;
            P[i * N + j] = 0.0;
            R[i * N + j] = 0.0;
        }
    }

    // Inicializamos la primer barrera donde tienen que pasar T elementos
    pthread_barrier_init(&barrera, NULL, T);

    // Inicio de la operación paralela
    double timetick = dwalltime();

    // Le damos comportamiento y creamos los hilos
    for (int i = 0; i < T; i++) {
        threads_ids[i] = i;
        pthread_create(&hilos[i], NULL, &behavior, (void*)&threads_ids[i]);
    }

    // Esperamos a que terminen los T hilos
    for (int i = 0; i < T; i++) {
        pthread_join(hilos[i], NULL);
    }

    pthread_barrier_destroy(&barrera);

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

    // Liberamos memoria
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