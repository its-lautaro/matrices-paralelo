#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void printMatrixBlock(int* matrixBlock, int numRows, int numCols) {
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            printf("%d ", matrixBlock[i * numCols + j]);
        }
        printf("\n");
    }
}

void producto_matricial(int* A, int* B, int* AB, int N, int BS) {
    for (int i = 0; i < BS; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                AB[i * N + j] = A[i * N + k] * B[j * N + k];
            }

        }

    }

}

int main(int argc, char** argv) {
    int id, T;
    int N = 16; // Matrix dimension
    int BS = 4; // Block size

    //Este programa solo funciona si todos los procesos procesan la misma cantidad de bloques

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &T);

    int elementos_por_fila = N * BS;
    int displs[T]; // Desplazamientos de buffer de cada proceso
    int sendcounts[elementos_por_fila]; //cantidad de elementos que se envian en cada scatterv (los que haya en una fila de bloques)
    int test[T];

    for (int i = 0; i < T; i++)
    {
        test[i]=0;
    }
    

    int* A, * B, * AB;
    int* Ablk = (int*)malloc(BS * N * sizeof(int));; // Block of the matrix for each process
    int* ABblk = (int*)malloc(BS * N * sizeof(int));; // Block of the matrix for each process

    B = (int*)malloc(N * N * sizeof(int)); //Todos los procesos necesitan la matriz B

    if (id == 0) {
        // Allocate and initialize the matrix
        A = (int*)malloc(N * N * sizeof(int));
        AB = (int*)malloc(N * N * sizeof(int));

        for (int i = 0; i < N * N; i++) {
            A[i] = i + 1;
            B[i] = 1;
        }

        for (int i = 0; i < T; i++) {
            sendcounts[i] = elementos_por_fila;
        }
        // Calcula la posicion del buffer desde la que comienza a enviar datos el emisor, tiene en cuenta la id del proceso y la cantidad de filas que se le mandaron
        for (int j = 0; j < T; j++) {
            displs[j] = j * (N * BS);
        }
    }
    //Todos los procesos necesitan la matriz B entera
    MPI_Bcast(B, N * N, MPI_INT, 0, MPI_COMM_WORLD);
    int filas_por_proceso = (N / (BS * T)); //cantidad de filas de bloques por proceso

    while (filas_por_proceso > 0) {
        filas_por_proceso--;
        // Recibo trabajo
        //int MPI_Scatterv(const void *sendbuf,const int *sendcounts,const int *displs,MPI_Datatype sendtype,void *recvbuf,int recvcount,MPI_Datatype recvtype,int rootRank,MPI_Comm comm)
        MPI_Scatterv(A, sendcounts, displs, MPI_INT, Ablk, BS * N, MPI_INT, 0, MPI_COMM_WORLD);
        // Hago trabajo
        //printf("Process %d received matrix block:\n", id);
        //printMatrixBlock(Ablk, BS, N);
        producto_matricial(Ablk, B, ABblk, N, BS);
        // El proceso root calcula el nuevo desplazamiento
        //int MPI_Gatherv(const void *sendbuf, const int *sendcounts,MPI_Datatype sendtype,void *recvbuf,int *recvcounts,const int *displs,MPI_Datatype recvtype,int rootRank,MPI_Comm comm)
        MPI_Gatherv(ABblk, N * BS, MPI_INT, AB, sendcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
        if (id == 0) {
            for (int j = 0; j < T; j++) {
                displs[j] += 2 * (N * BS);
            }
        }
    }
    if (id == 0) {
        printf("Matriz reconstruida\n");
        printMatrixBlock(AB, N, N);
    }

    MPI_Finalize();
    return 0;
}
