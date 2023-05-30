#include <stdio.h>
#include <stdlib.h>  
#include <sys/time.h>
#include <pthread.h>

//variables compartidas
int N;//tamaño de la matriz
int BS;//tamaño de bloque
int T;//cantidad de threads
double* A, * B, * C, * D; //matrices
double* P, * ABC, * DCB, * DC, * AB, * R; //matrices resultado
double* maximos, * minimos, * sumas;
double maxD, minA;
double promP;

//declaracion de barreras
pthread_barrier_t barrera;

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

void validarAB(double* matriz) {
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

void validarABC(double* matriz) {
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

void validarP(double* matriz) {
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

void validarR(double* matriz) {
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
    minimos[id] = local_min;
    maximos[id] = local_max;
    pthread_barrier_wait(&barrera); //esperar a que se terminen de escribir las variables compartidas
    //reducir arreglo
    int procesos_activos = T / 2;
    while (procesos_activos > 0) {
        if (id < procesos_activos) {
            local_min = (minimos[id * 2] > minimos[(id * 2) + 1]) ? minimos[(id * 2) + 1] : minimos[id * 2];
            local_max = (maximos[id * 2] < maximos[(id * 2) + 1]) ? maximos[(id * 2) + 1] : maximos[id * 2];
        }
        pthread_barrier_wait(&barrera); //esperan todos a que se termine de leer las variables compartidas
        if (id < procesos_activos) {
            //se reduce el arreglo
            maximos[id] = local_max;
            minimos[id] = local_min;
        }
        procesos_activos /= 2; //reduce a la mitad los procesos activos
        pthread_barrier_wait(&barrera); //esperar a que se terminen de escribir las variables compartidas
    }
    //pasar los datos a maxD y minA
    if (id == 0) {
        maxD = maximos[0];
        minA = minimos[0];
    }
}

void sum_promedio(int id) {
    double* Pblk, * ABCblk, * DCBblk, * Ablk, * Bblk, * Cblk, * Dblk, * ABblk, * DCblk;
    double sum_local;
    //for (int I = (id*(N/T)); I < ((id+1)*(N/T)); I ++) { otra opcion 
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
    //cada proceso guarda su suma local
    sumas[id] = sum_local;
    pthread_barrier_wait(&barrera); //esperar a que se terminen de escribir las variables compartidas
    //reducir arreglo
    int procesos_activos = T / 2;
    while (procesos_activos > 0) {
        if (id < procesos_activos) {
            sum_local = 0;
            sum_local = sumas[id * 2] + sumas[id * 2 + 1];
        }
        pthread_barrier_wait(&barrera);
        if (id < procesos_activos) {
            sumas[id] = sum_local;
        }
        procesos_activos /= 2;
        pthread_barrier_wait(&barrera);
    }
    if (id == 0) {
        promP = sumas[0] / (N * N);
    }
}

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

    //Calcula A.B.C y D.C.B, MinA y MaxD
    multBloques(id);

    pthread_barrier_wait(&barrera); //esperando MaxD y MinA

    // Calcula P = MaxD.ABC + MinA.DCB y PromP
    sum_promedio(id);

    pthread_barrier_wait(&barrera); //esperando promP

    //Calcula R = PromP.P
    producto_escalar(id);

    pthread_exit(NULL);
}

int main(int argc, char* argv[]) {
    BS=64; //tamaño de bloque optimo
    
    // Chequeo de parámetros
    if ((argc != 3) || ((N = atoi(argv[1])) <= 0) || ((T = atoi(argv[2])) <= 0)) {
        printf("Error en los parámetros. Usar: ./%s N T\n", argv[0]);
        exit(1);
    }

    printf("Matriz %dx%d en %d hilos\n",N,N,T);
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
            A[i * N + j] = (i * N + j) + 1;
            B[j * N + i] = 1.0; //ordenada por columnas
            AB[i * N + j] = 0.0;
            C[j * N + i] = 1.0; //ordenada por columnas
            D[i * N + j] = (i * N + j) + 1;
            DC[i * N + j] = 0.0;
            ABC[i * N + j] = 0.0;
            DCB[i * N + j] = 0.0;
            P[i * N + j] = 0.0;
            R[i * N + j] = 0.0;
        }
    }

    //Tomar tiempo start
    double timetick = dwalltime();

    //Inicializamos la primer barrera donde tienen que pasar T elementos
    pthread_barrier_init(&barrera, NULL, T);

    //Le damos comportamiento y creamos los hilos
    for (int i = 0; i < T; i++) {
        threads_ids[i] = i;
        pthread_create(&hilos[i], NULL, &behavior, (void*)&threads_ids[i]);
    }

    //Esperamos a que terminen los T hilos
    for (int i = 0; i < T; i++) {
        pthread_join(hilos[i], NULL);
    }

    pthread_barrier_destroy(&barrera);

    double totalTime = dwalltime() - timetick;
    printf("Tiempo en bloques de %d x %d: %f\n", BS, BS, totalTime);


    printf("Maximo:%f Minimo:%f Promedio:%f\n", maxD, minA, promP);
    //Se validan los resultados para una matriz no simetrica (cuyos elementos son (i*N+j))
    printf("Validando producto AB.. ");
    validarAB(AB);
    printf("Validando producto DC.. ");
    validarAB(DC);
    printf("Validando producto ABC.. ");
    validarABC(ABC);
    printf("Validando producto DCB.. ");
    validarABC(DCB);
    //maximo (N*N) minimo 1
    printf("Validando matriz P.. ");
    validarP(P);
    printf("Validando matriz R.. ");
    validarR(R);

    //Liberamos memoria
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