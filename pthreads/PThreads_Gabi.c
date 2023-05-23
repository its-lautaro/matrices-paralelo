//PThreads 
//Tengo que hacer MaxD, sumBloques y promP

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

//variables compartidas
int N, T;
double maxValueD = -1, sumTotal = 0;
double* D;
double* P;
pthread_mutex_t mutex_maxD, mutex_promP;

//calculamos el valor maximo de D de forma secuencial
/*double maxD_s(double *matriz, int n){
	double maxD = -1;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (matriz[i * n + j] > maxD) {
                maxD = matriz[i * n + j];
            }
        }
    }
}*/

//calculamos la suma de dos matrices
/*void sumBloques(double* A, double* B, double* C, int n, int bs) {
    //suma en bloques
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
                            cblk[i * n + j] = ablk[i * n + k] + bblk[i * n + k];
                        }
                    }
                }
            }
        }
    }
    //suma 
    // for (int i = 0; i < n; i++)
    // {
    //     for (int j = 0; j < n; j++)
    //     {
    //         C[i*n+j] += A[i*n+j] + B[i*n+j];
    //     }

    // }   
}
*/

//calculamos el promedio de P
/*double promP_s(double *matriz, int n){
	int sum = 0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sum += P[i * n + j];
        }
    }

    return (sum / (n*n));
}*/

void* maxD(void *arg){
	//guardo mi propio id
	int tid=*(int*)arg;

	int comienzo = tid * (N*N/T);
	int fin = comienzo + (N*N/T); 
	double max = D[comienzo];

	//recorro la matriz de acuerdo a mi porcion asignada
	for(int i = comienzo; i < fin; i++){
		if(D[i] > max){
			max = D[i];
		}
	}

	//comparo por exclusion mutua mi mayor con la variable compartida maxValueD
	pthread_mutex_lock(&mutex_maxD);
		
	if(max > maxValueD){
		maxValueD = max;
	}

	pthread_mutex_unlock(&mutex_maxD);

 	pthread_exit(NULL);
}

void* promP(void *arg){
	//guardo mi propio id
	int tid=*(int*)arg;

	int comienzo = tid * (N*N/T);
	int fin = comienzo + (N*N/T);
	double sum = 0;
	
	//recorro la matriz de acuerdo a mi porcion asignada
	for(int i = comienzo; i < fin; i++){
		sum += P[i];
	}

	pthread_mutex_unlock(&mutex_promP);
		sumTotal += sum;
	pthread_mutex_lock(&mutex_promP);

	pthread_exit(NULL);
}

int main(int argc, char* argv[]){
	N = atoi(argv[1]);
	T = atoi(argv[2]);
	pthread_t misThreads[T];
	int threads_ids[T];

	//Alocacion
	D = (double*)malloc(N * N * sizeof(double));
	P = (double*)malloc(N * N * sizeof(double));

	//Inicializacion
	for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            D[i * N + j] = i*i;
            P[i * N + j] = i*i;
        }
    }

	pthread_mutex_init(&mutex_maxD, NULL);
	pthread_mutex_init(&mutex_promP, NULL);

	//asigno los identificadores y creo los hilos
	for(int id=0;id<T;id++){
		threads_ids[id]=id;
		pthread_create(&misThreads[id],NULL,&promP,(void*)&threads_ids[id]);
	}

	//espero hasta que terminen todos los hilos
	for(int id=0;id<T;id++){
		pthread_join(misThreads[id],NULL);
	}

	pthread_mutex_destroy(&mutex_maxD);
	pthread_mutex_destroy(&mutex_promP);

	double promedioP = sumTotal/(N*N);

	free(D);
	free(P);

	return 0;
}
