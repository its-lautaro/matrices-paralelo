//Busca el tamaño de bloque optimo
#include<stdio.h>
#include<stdlib.h>  
#include <sys/time.h>

double dwalltime() {
	double sec;
	struct timeval tv;

	gettimeofday(&tv, NULL);
	sec = tv.tv_sec + tv.tv_usec / 1000000.0;
	return sec;
}

int main(int argc, char* argv[]) {
	double* A, * B, * C;
	int n= atoi(argv[1]);
	int bs = 16;
	//bloques
	double* ablk, * bblk, * cblk;
	//indices
	int I, J, K;
	int i, j, k;

	double timetick;

	// Alocar  
	A = (double*)malloc(n * n * sizeof(double));
	B = (double*)malloc(n * n * sizeof(double));
	C = (double*)malloc(n * n * sizeof(double));

	// Inicializacion
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			A[i * n + j] = 1.0;
			B[j * n + i] = 1.0; //ordenada por columnas
			C[i * n + j] = 0.0;
		}
	}

	while (bs <= n) {

		//tomar tiempo start
		timetick = dwalltime();

		for (I = 0; I < n; I += bs) {
			for (J = 0; J < n; J += bs) {
				cblk = &C[I * n + J];
				for (K = 0; K < n; K += bs) {
					ablk = &A[I * n + K];
					bblk = &B[J * n + K];
					for (i = 0; i < bs; i++) {
						for (j = 0; j < bs; j++) {
							for (k = 0; k < bs; k++) {
								cblk[i * n + j] += ablk[i * n + k] * bblk[j * n + k];
							}
						}
					}
				}
			}
		}

		//calcular tiempo de ejecución
		double totalTime = dwalltime() - timetick;
		printf("Tiempo en bloques de %d x %d: %f\n", bs, bs, totalTime);
		bs*=2;
	}


	free(A);
	free(B);
	free(C);

	return 0;
}



