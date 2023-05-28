#include <stdio.h>
#include <stdlib.h>  
#include <mpi.h>
#include <sys/time.h>

double dwalltime() {
    double sec;
    struct timeval tv;

    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

void printMatriz(double* matriz, int N) {
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%f ", matriz[i * N + j]);
        }
        printf("\n");
    }
}

double sum_promedio(double* ABC, double* DCB, double* P, int N, double minA, double maxD, int BS, int comienzo, int fin) {
    double* Pblk, * ABCblk, * DCBblk;
    double sum = 0;

    //for (int I = (id*(N/T)); I < ((id+1)*(N/T)); I ++) { otra opcion 
    for (int I = comienzo; I < fin; I++) {
        for (int J = 0; J < N; J += BS) {
            Pblk = &P[I * N + J];
            ABCblk = &ABC[I * N + J];
            DCBblk = &DCB[I * N + J];
            for (int i = 0; i < BS; i++) {
                for (int j = 0; j < BS; j++) {
                    Pblk[i * N + j] = maxD * ABCblk[i * N + j] + minA * DCBblk[i * N + j];
                    sum += Pblk[i * N + j];
                }
            }
        }
    }

    return sum;
}

void producto_escalar(double* P, double* R, int N, double promP, int BS,int comienzo, int fin) {
    double* Pblk, * Rblk;

    for (int I = comienzo; I < fin; I++){
        for (int J = 0; J < N; J += BS){
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

int main(int argc, char** argv) {
    const int BS = 2; 
    int N = 8;
    int T = 6;
    double* A, * B, * C, * D;
    double* P, * R, * AB, * ABC, * DC, * DCB;
    double* P_part, *ABC_part, *DCB_part, *R_part;
    double* sumas, *minimos, *maximos;
    double maxD = 1, minA = 1;
    double promP = 0, sum_local = 0.0;
    int id, nro_procesos;

    MPI_Init(&argc, &argv); //Inicializacion de ambiente
    MPI_Comm_rank(MPI_COMM_WORLD, &id); //identificador (rank) de cada proceso
    MPI_Comm_size(MPI_COMM_WORLD, &nro_procesos); //obtenemos el numero de procesos
 
    // Alocar las matrices resultados, solo lo realiza el root 
    if(id == 0){
        P = (double*)malloc(N * N * sizeof(double));
        ABC = (double*)malloc(N * N * sizeof(double));
        DCB = (double*)malloc(N * N * sizeof(double));
        R = (double*)malloc(N * N * sizeof(double));
    }

    //Cada uno asigna solo su parte
    P_part = (double*)malloc(BS * BS * sizeof(double));
    ABC_part = (double*)malloc(BS * BS * sizeof(double));
    DCB_part = (double*)malloc(BS * BS * sizeof(double));
    R_part = (double*)malloc(BS * BS * sizeof(double));

    // Inicializacion
    if(id == 0){
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                ABC[i * N + j] = 0.0;
                DCB[i * N + j] = 0.0;
                P[i * N + j] = 0.0;
                R[i * N + j] = 0.0;
            }
        }
    }

    //tomar tiempo start
    //double timetick = dwalltime();

    //Esta parte se hara anteriormente y nos dara como resultado ABC y BCD
    //ademas tendremos los minimos y maximos de las matrices A y D respectivamente
    //multBloques(A,B,AB,C,ABC,D,DC,DCB,N,BS,id,T,minimos,maximos);

    MPI_Barrier(MPI_COMM_WORLD); //primer barrera de sincronizacion

    //Distribuyo las matrices ABC y DCB entre los procesos
    MPI_Scatter(ABC,(BS*BS),MPI_DOUBLE,ABC_part,(BS*BS),MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Scatter(DCB,(BS*BS),MPI_DOUBLE,DCB_part,(BS*BS),MPI_DOUBLE,0,MPI_COMM_WORLD);

    //P = MaxD*ABC + MinA*DCB, PromP
    sum_local = sum_promedio(ABC_part,DCB_part,P_part,N,minA,maxD,BS,0,(N/BS));

    //Recopilo los resultados parciales en el proceso raiz (P)
    MPI_Gather(P_part,(BS*BS),MPI_DOUBLE,P,(BS*BS),MPI_DOUBLE,0,MPI_COMM_WORLD);

    //Realizo la reduccion de las sumas
    MPI_Reduce(&sum_local,&promP,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD); //segunda barrera de sincronizacion

    if(id == 0){
        promP = promP/(N*N); 
        printf("El promedio es igual a :%f\n", promP);
    }

    //deberia hacer un broadcast para que todos tengan el promedio calculado por el root? 

    MPI_Barrier(MPI_COMM_WORLD); //tercer barrera de sincronizacion   

    //Distribuyo la matriz P entre los procesos
    MPI_Scatter(P,(BS*BS),MPI_DOUBLE,P_part,(BS*BS),MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    //R=promP*P
    producto_escalar(P_part,R_part,N,promP,BS,0,(N/BS));

    //Recopilo los resultados parciales en el proceso raiz (R)
    MPI_Gather(R_part,(BS*BS),MPI_DOUBLE,R,(BS*BS),MPI_DOUBLE,0,MPI_COMM_WORLD);



    MPI_Barrier(MPI_COMM_WORLD); //barrera de sincronizacion auxiliar para imprimir
    printMatriz(R,N);




    //double totalTime = dwalltime() - timetick;
    //printf("Tiempo en bloques de %d x %d: %f\n", BS, BS, totalTime);

    //liberamos memoria
    if(id == 0){
        free(ABC);
        free(DCB);
        free(P);
        free(R);
    }
    free(ABC_part);
    free(DCB_part);
    free(P_part);
    free(R_part);

    MPI_Finalize(); //Finalizacion del ambiente
    return 0;
}