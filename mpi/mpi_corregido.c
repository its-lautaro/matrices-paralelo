#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void printMatrixBlock(long double* matrixBlock, int numRows, int numCols) {
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            printf("%Lf ", matrixBlock[i * numCols + j]);
        }
        printf("\n");
    }
}


void producto_matricial(int* A, int* B, int* AB, int* D, int* C, int* DC, int* ABC, int* DCB, int* max, int* min, int N, int BS, int nro_procesos) {
    int filas = N / nro_procesos;
    int local_min = N * N;
    int local_max = -1;

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
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                ABC[i * N + j] += AB[i * N + k] * C[j * N + k];
                DCB[i * N + j] += DC[i * N + k] * B[j * N + k];
            }
        }
    }
    MPI_Reduce(&local_max, max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_min, min, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
}

double sum_promedio(int* ABC, int* DCB, int* P, int N, int minA, int maxD, long int* sumTotal, int BS, int nro_procesos) {
    int filas = N / nro_procesos;
    long sum = 0;

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
    //Realizo la reduccion de las sumas
    //MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int rootRank, MPI_Comm comm)
    MPI_Reduce(&sum, sumTotal, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
}


void producto_escalar(int* P, long double* R, int N, long double promP, int BS, int nro_procesos) {
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
    int id, nro_procesos;
    int N = 4; // Matrix dimension
    int BS = 2; // Block size
    int maxD, minA;
    long sumTotal;
    long double promP;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &nro_procesos);

    int* A, * B, * AB, * D, * C, * DC, * DCB, * ABC, * P;
    long double* R;
    int* Ablk = (int*)malloc(((N / nro_procesos) * N) * sizeof(int)); // Block of the matrix for each process
    int* ABblk = (int*)calloc(((N / nro_procesos) * N), sizeof(int)); // Block of the matrix for each process
    int* Dblk = (int*)malloc(((N / nro_procesos) * N) * sizeof(int)); // Block of the matrix for each process
    int* DCblk = (int*)calloc(((N / nro_procesos) * N), sizeof(int)); // Block of the matrix for each process
    int* Cblk = (int*)malloc(((N / nro_procesos) * N) * sizeof(int)); // Block of the matrix for each process
    int* DCBblk = (int*)calloc(((N / nro_procesos) * N), sizeof(int)); // Block of the matrix for each process
    int* ABCblk = (int*)calloc(((N / nro_procesos) * N), sizeof(int)); // Block of the matrix for each process
    int* Pblk = (int*)malloc(((N / nro_procesos) * N) * sizeof(int));
    long double* Rblk = (long double*)malloc(((N / nro_procesos) * N) * sizeof(long double));

    B = (int*)malloc(N * N * sizeof(int)); //Todos los procesos necesitan la matriz B
    C = (int*)malloc(N * N * sizeof(int)); //Todos los procesos necesitan la matriz C

    if (id == 0) {
        // Allocate and initialize the matrix
        A = (int*)malloc(N * N * sizeof(int));
        AB = (int*)malloc(N * N * sizeof(int));
        D = (int*)malloc(N * N * sizeof(int));
        DC = (int*)malloc(N * N * sizeof(int));
        P = (int*)malloc(N * N * sizeof(int));
        R = (long double*)malloc(N * N * sizeof(long double));

        for (int i = 0; i < N * N; i++) {
            A[i] = i + 1;
            B[i] = 1;
            D[i] = i + 1;
            C[i] = 1;
        }
    }

    //Todos los procesos necesitan la matriz B entera
    MPI_Bcast(B, N * N, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(C, N * N, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Scatter(A, ((N / nro_procesos) * N), MPI_INT, Ablk, ((N / nro_procesos) * N), MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(D, ((N / nro_procesos) * N), MPI_INT, Dblk, ((N / nro_procesos) * N), MPI_INT, 0, MPI_COMM_WORLD);

    producto_matricial(Ablk, B, ABblk, Dblk, C, DCblk, ABCblk, DCBblk, &maxD, &minA, N, BS, nro_procesos);

    MPI_Bcast(&minA, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&maxD, 1, MPI_INT, 0, MPI_COMM_WORLD);

    sum_promedio(ABCblk, DCBblk, Pblk, N, minA, maxD, &sumTotal, BS, nro_procesos);
    //Root calcula el promedio
    if (id == 0) promP = sumTotal / (N * N);
    MPI_Bcast(&promP, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);

    //R=promP*P
    producto_escalar(Pblk,Rblk,N,promP,BS,nro_procesos);

    //Recopilo los resultados parciales en el proceso raiz (R)
    MPI_Gather(Rblk,((N/nro_procesos)*N),MPI_LONG_DOUBLE,R,((N/nro_procesos)*N),MPI_LONG_DOUBLE,0,MPI_COMM_WORLD);

    if (id == 0) {
        printMatrixBlock(R,N,N);
    }

    MPI_Finalize();
    return 0;
}