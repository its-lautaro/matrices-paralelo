/* Encontrar el valor minimo en matrices de nxn
* Memoria Compartida (pthread)
* Lautaro
*/
#include<stdio.h>
#include<stdlib.h>  
#include <sys/time.h>
#include <pthread.h>

//variables compartidas
int N;//tamaño de la matriz
int T;//cantidad de threads
double* matriz;
double minimo=9999;
pthread_mutex_t mutex_minimo;

double dwalltime() {
    double sec;
    struct timeval tv;

    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

double sec_min(double* matriz, int n) {
    double min = 9999;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (matriz[i * n + j] < min) {
                min = matriz[i * n + j];
            }
        }
    }

    return min;
}

void* min(void* arg) {
    int id = *(int*)arg;
    int start = id * ((N*N) / T);
    int end = id < (T-1)? (id + 1) * ((N*N) / T):(N*N)-1;
    double min_local = matriz[start];

    for (int i = start + 1; i < end; i++) {
        if (matriz[i] < min_local) {
            min_local = matriz[i];
        }
    }

    pthread_mutex_lock(&mutex_minimo);
        if (min_local < minimo) {
            minimo = min_local;
        }
    pthread_mutex_unlock(&mutex_minimo);

    pthread_exit(NULL);
}

int main(int argc, char* argv[]) {
    N = atoi(argv[1]); //tamaño matriz
    T = atoi(argv[2]); //cantidad de threads

    double timetick;
    matriz = (double*)malloc(sizeof(double) * N * N);

    pthread_t hilos[T];
    int threads_ids[T];

    //inicializar
    pthread_mutex_init(&mutex_minimo, NULL);
    for (int i = 0;i <= (N * N);i++) {
        matriz[i] = i+1;
    }

    //tomar tiempo start
    timetick = dwalltime();

    for (int i = 0; i < T; i++) {
        threads_ids[i] = i;
        pthread_create(&hilos[i], NULL, &min, (void*)&threads_ids[i]);
    }
    for (int i = 0; i < T; i++) {
        pthread_join(hilos[i], NULL);
    }
    pthread_mutex_destroy(&mutex_minimo);

    double totalTime = dwalltime() - timetick;

    printf("Tiempo %fs\n", totalTime);
    printf("valor minimo: %f\n", minimo);
    return 0;
}