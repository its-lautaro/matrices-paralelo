# Sistemas Distribuidos y Paralelos: Trabajo de promoción

## Calculo del tamaño óptimo del bloque de memoria
Obtenemos el promedio de 3 ejecuciones utilizando distintos tamaños de bloque y distintos tamaños de matriz.

| N / BS   | 32      | 64      | 128         | 256         | 512        | 1024    | 2048    | 4096 |
| -------- | ------- | ------- | ----------- | ----------- | ---------- | ------- | ------- | ---- |
| **1024** | 13.207  | 13.038  | 12.873      | 12.797      | **12.761** | -       | -       | -    |
| **2048** | 105.212 | 104.202 | 102.841     | **102.329** | 105.564    | 104.719 | -       | -    |
| **4096** | 839.167 | 833.465 | **823.959** | 855.641     | 846.233    | 838.107 | 833.543 | -    |

Los tamaños de bloque que obtienen el mejor tiempo de ejecución difieren entre cada matriz, pero si comparamos la diferencia de tiempos de las matrices al usar un bloque de tamaño de 128 con respecto a su bloque óptimo la diferencia es prácticamente despreciable. Por este motivo elegimos el tamaño óptimo para todas las matrices en **128x128**.

>Promedio de 3 ejecuciones

#### Script utilizado para correr el programa
```bash
#!/bin/bash
#SBATCH --job-name=multBloques  # Job name
#SBATCH --array=32,64,128,256,512         # Parameter values to be used
#SBATCH --ntasks=1              # Number of tasks (or processes) per job array task

#SBATCH -N 1
#SBATCH --exclusive

# Define the parameter values
params=($SLURM_ARRAY_TASK_ID)

# Run the program with the current parameter value
output_file="/nethome/sdyp4/ejercicio1/1024/output${params[0]}.txt"  # Output file name based on the parameter value
/nethome/sdyp4/ejercicio1/multBloques 1024 ${params[0]} > $output_file
```
## Algoritmo secuencial

Cuando existe dependencia de datos entre las operaciones realizadas, resulta en mejores tiempos de ejecución realizar estas operaciones simultaneamente, ya que minimizamos los accesos a memoria y ademas aprovechamos que los datos ya se encuentran en cache. Las mejoras se pueden observar ejecutando los programas expMatrices1, 2 y 3 de la practica 1.

### expMatrices1.c: AB+AC+AD
El programa expMatrices1 realiza la operacion de dos formas:

* Suma los elementos de las matrices AB, AC yAD a medida que los calcula.
* Calcula las matrices completas, cada una en su propio bucle, luego las suma.

| N    | 1 bucle | 4 bucles |
| ---- | ------- | -------- |
| 1024 | 36.582s | 38.912s  |

Al haber dependencia de datos, realizar las operaciones en 1 bloque evita fallos de cache (se suman valores que ya estan en memoria)

### expMatrices2.c: AB+CD
Observamos lo mismo que en el programa anterior.

| N    | 1 bucle | 3 bucles |
| ---- | ------- | -------- |
| 1024 | 24.863s | 25.877s  |

### expMatrices3.c: BA+CAD
Observamos lo mismo que en el programa anterior.

| N    | 2 bucles | 4 bucles |
| ---- | -------- | -------- |
| 1024 | 37.577s  | 39.985s  |

### Probando alternativas
Teniendo en cuenta los resultados obtenidos, probaremos que esta tecnica es aplicable a nuestra operación, la misma consiste en los siguientes calculos

- R = PromP.P
- P = MaxD.(ABC) + MinA.(DCB)

Es decir, 4 multiplicaciones matriciales (AB,ABC,DC,DCB), cálculo de máximo y minimo, 3 multiplicaciones escalares (MaxD.ABC, MinA.DCB,PromP.P), una suma y el cálculo de un promedio. Una forma de aprovechar la dependencia de datos para mejorar los fallos de cache, podría ser:
- Realizar las multiplicaciones matriciales en su propio bucle, y al mismo tiempo calcular el máximo y mínimo de las matrices D y A respectivamente.
- Una vez calculados MaxD y MinD se puede realizar el calculo de P en otro bucle. Mientras llevamos acabo este cálculo dentro del mismo bloque podemos llevar a cabo el calculo de PromP (el promedio de P)
- Con el promedio calculado finalmente podemos realizar el cálculo de R.

Desarrollamos esta y otra alternativa que realiza cada operacion en su propio bucle para comparar los tiempos de ejecución. El tamaño de bloque utilizado fue el tamaño óptimo obtenido para cada N, en el punto 1.

**Realizando una operación por bucle**

| 1024   | 2048    | 4096     |
| ------ | ------- | -------- |
| 51.278 | 409.614 | 3294.489 |

**Realizando múltiples operaciones por bucle**
  
| 1024   | 2048    | 4096     |
| ------ | ------- | -------- |
| 50.226 | 408.717 | 3289.593 |

Comprobamos que aprovechar la dependencia de datos mejora el tiempo de ejecución. Utilizaremos el programa secuencial que realiza multiples operaciones por bucle para obtener la solución paralela.

## Algoritmo paralelo (Memoria compartidad)
Para paralelizar nuestro algoritmo secuencial, comenzamos distinguiendo el trabajo que debe llevarse a cabo. Como se describió en el algoritmo secuencial, nuestro trabajo puede dividirse en 3 tareas que conllevan:

- La multiplicacion matricial (ABC y DCB) y el cálculo de máximo y mínimo (min(A) y max(D))
- La multiplicacion escalar y la suma de matrices (P=MinA.ABC + MaxD.DCB), y el cálculo del promedio (suma de los elementos P y division por cantidad de elementos en P)
- Multiplicación escalar (PromP.P)

Estas operaciones deben realizarse de forma secuencial, pero el trabajo de cada una puede paralelizarse ya que trabajamos con las matrices en bloques.

Ahora si analizamos la comunicación de los procesos, resulta que solo hay una instancia de comunicacion y es cuando se calculan el maximo y el minimo. Para este problema, que es un problema de reduccion, utilizaremos barreras.


### Pthread

### OpenMP

## Algoritmo paralelo (Memoria distribuida)

### MPI