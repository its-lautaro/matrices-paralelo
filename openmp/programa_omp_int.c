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

void validarAB(int* matriz, int N) {
    for (int i = 0; i < N; i++) {
        int res = 0; //valor utilizado para validar
        for (int j = 1; j < N + 1; j++) {
            res += i * N + j;
        }
        for (int j = 0; j < N; j++) {
            if (matriz[i * N + j] != res) {
                printf("ERROR %d no coincide con %d\n", matriz[i * N + j], res);
                return;
            }
        }
    }
    printf("OK\n");
}

void validarABC(int* matriz, int N) {
    int res;
    for (int i = 0; i < N; i++) {
        res = 0; //valor utilizado para validar
        for (int j = 0; j < N; j++) {
            res += (i * N + j) + 1; //elemento de AB fila i
        }
        res *= N; // elemento de ABC fila i
        for (int j = 0; j < N; j++) {
            if (matriz[i * N + j] != res) {
                printf("ERROR %d no coincide con %d\n", matriz[i * N + j], res);
                return;
            }
        }
    }
    printf("OK\n");
}

void validarP(int* matriz, int N) {
    int res;
    for (int i = 0; i < N; i++) {
        res = 0; //valor utilizado para validar
        for (int j = 0; j < N; j++) {
            res += (i * N + j) + 1; //elemento de AB fila i
        }
        res *= N; // elemento de ABC fila i
        res = res * 1 + res * (N * N); //elemento de P fila i
        for (int j = 0; j < N; j++) {
            if (matriz[i * N + j] != res) {
                printf("ERROR %d no coincide con %d\n", matriz[i * N + j], res);
                return;
            }
        }
    }
    printf("OK\n");
}

void validarR(int* matriz, int N, long double promP) {
    int res;
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
                printf("ERROR %d no coincide con %d\n", matriz[i * N + j], res);
                return;
            }
        }
    }
    printf("OK\n");
}

int main() {
    int N = 4096; // Suponiendo matrices cuadradas de tamaño N x N
    int BS = 64; // Tamaño del bloque
    int T = 6; // Número de hilos a utilizar

    int* A, * B, * AB, * C, * ABC, * D, * DCB, * DC, * P, * R;
    int maximos[T];
    int minimos[T];
    int maxD = -1, minA = 9999;
    
    long sumas[T];
    long sumaTotal = 0;
    long double promP;

    // Asignar memoria para las matrices A, B y AB
    A = (int*)malloc(N * N * sizeof(int));
    B = (int*)malloc(N * N * sizeof(int));
    C = (int*)malloc(N * N * sizeof(int));
    AB = (int*)malloc(N * N * sizeof(int));
    DC = (int*)malloc(N * N * sizeof(int));
    ABC = (int*)malloc(N * N * sizeof(int));
    D = (int*)malloc(N * N * sizeof(int));
    DCB = (int*)malloc(N * N * sizeof(int));
    P = (int*)malloc(N * N * sizeof(int));
    R = (int*)malloc(N * N * sizeof(int));

    // Inicializar las matrices A, B y AB
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i * N + j] = (i * N + j) + 1;
            B[j * N + i] = 1;
            AB[i * N + j] = 0;
            C[j * N + i] = 1;
            ABC[i * N + j] = 0;
            D[i * N + j] = (i * N + j) + 1;
            DC[i * N + j] = 0;
            DCB[i * N + j] = 0;
        }
    }
    double timetick = dwalltime();
    // Realizar la multiplicación de matrices en bloques (ABC y DCB) y construir arreglo de maximos y minimos
#pragma omp parallel num_threads(T)
    {
        int id = omp_get_thread_num();
        int min_local = 9999, max_local = -1;
#pragma omp for
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
                            //MIN(A)
                            if (A[i * N + j] < min_local) {
                                min_local = A[i * N + j];
                            }

                            //MAX(D)
                            if (D[i * N + j] > max_local) {
                                max_local = D[i * N + j];
                            }
                        }
                    }
                }
            }
            //Una vez que calcule un bloque de AB, calculo un bloque de ABC
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
        minimos[id] = min_local;
        maximos[id] = max_local;
    }
    //Barrera para esperar que todos los procesos terminen de construir los arreglos maximos y minimos
#pragma omp parallel for reduction(max: maxD) reduction(min: minA)
    for (int i = 0; i < T; i++) {
        if (maximos[i] > maxD) maxD = maximos[i];
        if (minimos[i] < minA) minA = minimos[i];
    }
    //Calcular la suma P=maxD*ABC + minA*DCB, en bloques y crear el arreglo sumas
#pragma omp parallel num_threads(T)
    {
        long suma_local = 0;
        int id = omp_get_thread_num();
#pragma omp for
        for (int I = 0; I < N; I += BS) {
            for (int J = 0; J < N; J += BS) {
                // Calcular la suma en bloques
                for (int i = I; i < I + BS; i++) {
                    for (int j = J; j < J + BS; j++) {
                        P[i * N + j] = maxD * ABC[i * N + j] + minA * DCB[i * N + j];
                        suma_local += P[i * N + j];
                    }
                }
            }
        }

        //printf("%d: suma local: %ld\n", id, suma_local);
        sumas[id] = suma_local;
    }
    //Barrera para esperar que todos los procesos terminen de construir el arreglo sumas
#pragma omp parallel for reduction(+: sumaTotal)
    for (int i = 0; i < T; i++) {
        sumaTotal += sumas[i];
    }
    //printf("Suma total %ld\n",sumaTotal);
    //El proceso maestro calcula promP
    promP = sumaTotal / (N * N);
    //Calcula el producto escalar R=PromP*P
#pragma omp parallel num_threads(T)
    {
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

    double totalTime = dwalltime() - timetick;
    printf("Tiempo en bloques de %d x %d: %f\n", BS, BS, totalTime);

    printf("Maximo:%d Minimo:%d Promedio:%Lf\n", maxD, minA, promP);
    //Se validan los resultados para una matriz no simetrica (cuyos elementos son (i*N+j))
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
