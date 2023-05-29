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

double sum_promedio(double* ABC, double* DCB, double* P, int N, double minA, double maxD, int BS, int id, int nro_procesos) {
    int filas = N/nro_procesos;
    double sum = 0;

    for (int I = 0; I < filas; I += BS) {
        for (int J = 0; J < N; J += BS) {
            printf("id:%d, I: %d, J:%d.\n", id, I,J);
            for (int i = I; i < I + BS; i++) {
                for (int j = J; j < J + BS; j++) {
                    P[i * N + j] = maxD * ABC[i * N + j] + minA * DCB[i * N + j];
                    printf("id:%d, ABC+BCD=%f+%f=P=%f.\n", id, ABC[i * N + j],DCB[i * N + j],P[i * N + j]);
                    sum += P[i * N + j];
                }
            }
        }
    }

    return sum;
}

void producto_escalar(double* P, double* R, int N, double promP, int BS, int nro_procesos) {
    int filas = N/nro_procesos;
    for (int I = 0; I < filas; I += BS){
        for (int J = 0; J < N; J += BS){
            for (int i = I; i < I + BS; i++) {
                for (int j = J; j < J + BS; j++) {
                    R[i * N + j] = promP * P[i * N + j];
                }
            }
        }
    }
}

int main(int argc, char** argv) {
    const int BS = 2; 
    int N = 8;
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
    P_part = (double*)malloc(((N/nro_procesos)*N) * sizeof(double));
    ABC_part = (double*)malloc(((N/nro_procesos)*N) * sizeof(double));
    DCB_part = (double*)malloc(((N/nro_procesos)*N) * sizeof(double));
    R_part = (double*)malloc(((N/nro_procesos)*N) * sizeof(double));

    // Inicializacion
    if(id == 0){
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                ABC[i * N + j] = i*N+j;
                DCB[i * N + j] = i*N+j;
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

    //MPI_Barrier(MPI_COMM_WORLD); //primer barrera de sincronizacion

    //Distribuyo las matrices ABC y DCB entre los procesos
    //MPI_Scatter(void *sendbuf,   int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int rootRank, MPI_Comm comm)
    MPI_Scatter(ABC,((N/nro_procesos)*N),MPI_DOUBLE,ABC_part,((N/nro_procesos)*N),MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Scatter(DCB,((N/nro_procesos)*N),MPI_DOUBLE,DCB_part,((N/nro_procesos)*N),MPI_DOUBLE,0,MPI_COMM_WORLD);

    //P = MaxD*ABC + MinA*DCB, PromP
    sum_local = sum_promedio(ABC_part,DCB_part,P_part,N,minA,maxD,BS,id,nro_procesos);

    printf("id: %d sume :%f\n", id, sum_local);

    //Recopilo los resultados parciales en el proceso raiz (P)
    MPI_Gather(P_part,((N/nro_procesos)*N),MPI_DOUBLE,P,((N/nro_procesos)*N),MPI_DOUBLE,0,MPI_COMM_WORLD);

    //MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int rootRank, MPI_Comm comm)
    //Realizo la reduccion de las sumas
    MPI_Reduce(&sum_local,&promP,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD); //segunda barrera de sincronizacion

    if(id == 0){
        promP = promP/(N*N); 
    }
    
    //MPI_Bcast(void *buffer, int count, MPI_Datatype dtype, int root, MPI_Comm comm)
    MPI_Bcast(&promP, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    printf("SOY: %d El promedio es igual a :%f\n", id, promP);

    MPI_Barrier(MPI_COMM_WORLD); //tercer barrera de sincronizacion   

    //Distribuyo la matriz P entre los procesos
    MPI_Scatter(P,((N/nro_procesos)*N),MPI_DOUBLE,P_part,((N/nro_procesos)*N),MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    //R=promP*P
    producto_escalar(P_part,R_part,N,promP,BS,nro_procesos);

    //Recopilo los resultados parciales en el proceso raiz (R)
    MPI_Gather(R_part,((N/nro_procesos)*N),MPI_DOUBLE,R,((N/nro_procesos)*N),MPI_DOUBLE,0,MPI_COMM_WORLD);

    double totalTime = dwalltime() - timetick;
    printf("Tiempo en bloques de %d x %d: %f\n", BS, BS, totalTime);

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