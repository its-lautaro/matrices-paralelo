### SISTEMAS DISTRIBUIDOS Y PARALELOS


## **PARALELIZACIÓN DE ALGORITMOS**

# TRABAJO DE PROMOCIÓN 2023

### Grupo 4
### **Integrantes**

* La Vecchia, Lautaro - 02031/2
* Torres, Gabriel - 1987/4
<br/>

## Indice
- [Introducción](#introducción)
- [Desarrollo](#desarrollo)
  - [Problema 1: C = AB, siendo A, B y C matrices cuadradas de N\*N](#problema-1-c--ab-siendo-a-b-y-c-matrices-cuadradas-de-nn)
    - [Cálculo del tamaño óptimo del bloque de memoria](#cálculo-del-tamaño-óptimo-del-bloque-de-memoria)
  - [Problema 2: R = PromP \* P, con P = MaxD \* (A*B*C) + MinA \* (D*C*B)](#problema-2-r--promp--p-con-p--maxd--abc--mina--dcb)
    - [Solución secuencial.](#solución-secuencial)
      - [expMatrices1.c: AB+AC+AD](#expmatrices1c-abacad)
      - [expMatrices2.c: AB+CD](#expmatrices2c-abcd)
      - [expMatrices3.c: BA+CAD](#expmatrices3c-bacad)
    - [Estrategias de implementación.](#estrategias-de-implementación)
    - [](#)
    - [Paralelización del algoritmo secuencial.](#paralelización-del-algoritmo-secuencial)
      - [PTHREADS](#pthreads)
      - [OPENMP](#openmp)
      - [MPI](#mpi)
- [Análisis de escalabilidad y conclusiones.](#análisis-de-escalabilidad-y-conclusiones)
    - [Tiempo de ejecución: Pthreads](#tiempo-de-ejecución-pthreads)
    - [Tiempo de ejecución: OpenMP](#tiempo-de-ejecución-openmp)
    - [Tiempo de ejecución: MPI](#tiempo-de-ejecución-mpi)
  - [SPEEDUP](#speedup)
    - [SPEEDUP para Pthreads](#speedup-para-pthreads)
    - [SPEEDUP para OpenMP](#speedup-para-openmp)
    - [SPEEDUP para MPI](#speedup-para-mpi)
  - [EFICIENCIA](#eficiencia)
    - [EFICIENCIA para Pthreads](#eficiencia-para-pthreads)
    - [EFICIENCIA para OpenMP](#eficiencia-para-openmp)
    - [EFICIENCIA para MPI](#eficiencia-para-mpi)
  - [ESCALABILIDAD](#escalabilidad)
    - [Observando los resultados obtenidos para Pthreads y OpenMP](#observando-los-resultados-obtenidos-para-pthreads-y-openmp)
    - [Observando los resultados obtenidos para MPI](#observando-los-resultados-obtenidos-para-mpi)

## Introducción

En el presente informe se explica la estrategia de diseño de los algoritmos que resuelven las operaciones con matrices y su posterior paralelizacion. 

Abordaremos por separado tres de las librerías principales para el manejo de paralelismo: **<span style="text-decoration:underline;">PThreads</span>**, **<span style="text-decoration:underline;">OpenMP</span>** y **<span style="text-decoration:underline;">MPI</span>**.  

El objetivo buscado es mejorar los tiempos de ejecución del algoritmo secuencial y posteriormente realizar un análisis de escalabilidad y eficiencia sobre las soluciones paralelas.

## Desarrollo

### Problema 1: C = AB, siendo A, B y C matrices cuadradas de N*N

Se implementó la multiplicación de dos matrices cuadradas utilizando la [técnica por bloques](#bookmark=id.nw952cvw80k3) guardando el resultado en otra matriz. El objetivo es obtener el tamaño óptimo de bloque para los distintos tamaños de matrices. 

El procedimiento seguido para obtener el tamaño óptimo consistió en probar distintos tamaños de bloque para distintos tamaños de matriz y medir el tiempo de ejecución, eligiendo aquellos que fueran más rápidos. A modo de ejemplo se ilustra uno de los scripts utilizados que automatizan la prueba de distintos tamaños de bloque para una matriz de 1024 elementos.


```
#!/bin/bash
#SBATCH --job-name=bs_optimo  # Nombre del job
#SBATCH --array=1-5             # Cantidad de jobs a generar
#SBATCH --ntasks=1              # Cantidad de procesos por job
#SBATCH -N 1
#SBATCH --exclusive

# Definir el arreglo de tamaños de bloque
bs=(32 64 128 256 512)

# Calcula el tamaño de bloque según el id del job
param_index=$((SLURM_ARRAY_TASK_ID-1))
bs_value=${bs[$param_index]}

# Ejecutar el programa con el valor actual de tamaño de bloque y guardar el output
output_file="1024/bs${bs_value}.txt"  # Nombre del archivo de salida
./mult_bloques 1024 ${bs_value} > $output_file
```

_<span style="text-decoration:underline;">Técnica por bloques:</span> estrategia de programación que consiste en dividir una tarea o proceso en bloques más pequeños y realizar las tarea en cada bloque de forma independiente. Nos permite mejorar el rendimiento al reducir el número de accesos a memoria y aprovechar la localidad de los datos._

#### Cálculo del tamaño óptimo del bloque de memoria

Obtenemos el promedio del tiempo de ejecución a partir de 3 ejecuciones independientes, utilizando distintos tamaños de bloque y distintos tamaños de matriz.

<table>
  <tr>
   <td><strong>N/BS</strong>
   </td>
   <td><strong>32</strong>
   </td>
   <td><strong>64</strong>
   </td>
   <td><strong>128</strong>
   </td>
   <td><strong>256</strong>
   </td>
   <td><strong>512</strong>
   </td>
   <td><strong>1024</strong>
   </td>
   <td><strong>2048</strong>
   </td>
   <td><strong>4096</strong>
   </td>
  </tr>
  <tr>
   <td><strong>1024</strong>
   </td>
   <td>13.207
   </td>
   <td>13.038
   </td>
   <td>12.873
   </td>
   <td>12.797
   </td>
   <td><strong>12.761</strong>
   </td>
   <td>-
   </td>
   <td>-
   </td>
   <td>-
   </td>
  </tr>
  <tr>
   <td><strong>2048</strong>
   </td>
   <td>105.212
   </td>
   <td>104.202
   </td>
   <td>102.841
   </td>
   <td><strong>102.329</strong>
   </td>
   <td>105.564
   </td>
   <td>104.719
   </td>
   <td>-
   </td>
   <td>-
   </td>
  </tr>
  <tr>
   <td><strong>4096</strong>
   </td>
   <td>839.167
   </td>
   <td>833.465
   </td>
   <td><strong>823.959</strong>
   </td>
   <td>855.959
   </td>
   <td>855.641
   </td>
   <td>846.233
   </td>
   <td>838.107
   </td>
   <td>-
   </td>
  </tr>
</table>

Los tamaños de bloque que obtienen el mejor tiempo de ejecución difieren entre cada matriz, pero si comparamos la diferencia de tiempos de las matrices al usar un bloque de tamaño de 128 (en fondo verde) con respecto a su bloque óptimo (en letra azul), la diferencia es prácticamente despreciable. Por este motivo elegimos el tamaño de bloque de 128 como el tamaño óptimo para todas las matrices, _siempre que sea posible*_.

> *Dependiendo de la cantidad de elementos que se reparta a cada proceso, puede darse que el tamaño de bloque máximo se vea limitado. En particular, al trabajar con 8 y 16 procesos se tiene que los tamaños de bloque serán:

<table>
  <tr>
   <td><strong>Unid_Proc/N</strong>
   </td>
   <td><strong>1024</strong>
   </td>
   <td><strong>2048</strong>
   </td>
   <td><strong>4096</strong>
   </td>
  </tr>
  <tr>
   <td><strong>8</strong>
   </td>
   <td><strong>64</strong>
   </td>
   <td>128
   </td>
   <td>128
   </td>
  </tr>
  <tr>
   <td><strong>16</strong>
   </td>
   <td><strong>32</strong>
   </td>
   <td><strong>64</strong>
   </td>
   <td>128
   </td>
  </tr>
</table>

### Problema 2: R = PromP * P, con P = MaxD * (A*B*C) + MinA * (D*C*B)

Cuando existe dependencia de datos entre las operaciones realizadas, como en este caso, resulta en mejores tiempos de ejecución efectuar estas operaciones simultáneamente, ya que se minimizan los accesos a memoria y además se aprovecha que los datos ya se encuentran en memoria caché. Las mejoras se pueden observar analizando el tiempo de ejecución de los programas expMatrices 1, 2 y 3 de la práctica 1.

#### Solución secuencial.
##### expMatrices1.c: AB+AC+AD

El programa expMatrices1 realiza la operación de dos formas:

* Suma los elementos de las matrices **AB**, **AC** y **AD** a medida que los calcula.
* Calcula las matrices **AB**, **AC** y **AD **completas , cada una en su propio bucle, y luego en otro bucle las suma.    

<table>
  <tr>
   <td>
<strong>N</strong>
   </td>
   <td><strong>1 bucle </strong>
   </td>
   <td><strong>4 bucles</strong>
   </td>
  </tr>
  <tr>
   <td>1024
   </td>
   <td>36.582s
   </td>
   <td>38.912s
   </td>
  </tr>
</table>

Al haber dependencia de datos, realizar las operaciones en 1 bucle evita fallos de cache (se suman valores que ya están en memoria)

##### expMatrices2.c: AB+CD

Observamos lo mismo que en el programa anterior.

<table>
  <tr>
   <td><strong>N</strong>
   </td>
   <td><strong>1 bucle </strong>
   </td>
   <td><strong>3 bucles</strong>
   </td>
  </tr>
  <tr>
   <td>1024
   </td>
   <td>24.863s
   </td>
   <td>25.877s
   </td>
  </tr>
</table>

##### expMatrices3.c: BA+CAD

Observamos lo mismo que en el programa anterior.

<table>
  <tr>
   <td><strong>N</strong>
   </td>
   <td><strong>2 bucles </strong>
   </td>
   <td><strong>4 bucles</strong>
   </td>
  </tr>
  <tr>
   <td>1024
   </td>
   <td>37.577s
   </td>
   <td>39.985s
   </td>
  </tr>
</table>

#### Estrategias de implementación.

Teniendo en cuenta los resultados obtenidos y aprovechando la técnica por bloques probaremos que esta técnica es aplicable a nuestra operación, la cual consiste en los siguientes cálculos:

* P = MaxD.(ABC) + MinA.(DCB)
* R = PromP.P

Es decir, 4 multiplicaciones matriciales (A\*B,AB\*C,D\*C,DC\*B), cálculo de máximo y mínimo, 3 multiplicaciones escalares (MaxD\*ABC, MinA\*DCB,PromP\*P), una suma y el cálculo de un promedio. Una forma de aprovechar la dependencia de datos para minimizar los fallos de caché, es dividir el cálculo en la menor cantidad de etapas posible:

1. _multBloques(): _Realizar las multiplicaciones matriciales en un mismo bucle, y al mismo tiempo calcular el máximo y mínimo de las matrices D y A respectivamente.
2. _sumPromedio()_: Una vez calculados MaxD y MinD se puede realizar el cálculo de P en otro bucle. Mientras se realiza este cálculo, dentro del mismo bloque se puede llevar a cabo la suma de todos los elementos de P, para el posterior cálculo de su promedio (PromP).
3. _productoEscalar(): _Con el promedio calculado finalmente podemos realizar el cálculo de R (P.PromP).

Se implementaron dos algoritmos secuenciales: 

* _secuencial_dependencia_datos.c_
* _secuencial.c_

 El primero utiliza la técnica anteriormente descrita. El segundo utiliza otra técnica alternativa en donde cada operación se realiza en su propio bucle. Comparando ambos tiempos de ejecución se comprobó que el tiempo es menor al realizar varias operaciones por bucle. El tamaño de bloque utilizado fue el tamaño óptimo (**[BS = 128](#cálculo-del-tamaño-óptimo-del-bloque-de-memoria-5)**).

<table>
  <tr>
   <td><strong>Bucles/N</strong>
   </td>
   <td><strong>1024</strong>
   </td>
   <td><strong>2048</strong>
   </td>
   <td><strong>4096</strong>
   </td>
  </tr>
  <tr>
   <td><strong>1 Op x Bucle</strong>
   </td>
   <td>51.641s
   </td>
   <td>409.614s
   </td>
   <td>3294.489
   </td>
  </tr>
  <tr>
   <td><strong>Múltiples Op x Bucle</strong>
   </td>
   <td>49.532s
   </td>
   <td>408.717s
   </td>
   <td>3289.593s
   </td>
  </tr>
</table>

Comprobamos que aprovechar la dependencia de datos mejora el tiempo de ejecución. Posteriormente se utilizará este programa secuencial (_secuencial_dependencia_datos.c_) para obtener la solución paralela.

#### Paralelización del algoritmo secuencial.

Con el objetivo de aprovechar al máximo las ventajas de utilizar técnica por bloques, se utiliza un enfoque para la paralelización en donde las operaciones internas siguen funcionando de forma secuencial, pero se paraleliza el bucle externo. Para lograr nuestro cometido, se utilizan las librerías anteriormente mencionadas optimizando sus diversas funcionalidades y estrategias. 

##### PTHREADS

Para implementar el comportamiento de los T procesos, se utilizaron barreras que sincronizan y separan las distintas etapas anteriormente descritas. Esto permite maximizar la eficiencia en términos de rendimiento y escalabilidad, ya que el sistema ahorra el costo en recursos asociado a crear y destruir hilos constantemente para cada etapa.

A su vez, dentro de cada etapa, se utilizó el id de cada proceso en el bucle externo para paralelizar las tareas, de esta manera cada proceso se reparte una fila de bloques de cada matriz para realizar las operaciones. Además, aprovechando que se recorren las matrices al realizar las operaciones, los procesos construyen arreglos para el cálculo de _MinA, MaxD y SumTotal (la suma total de elementos de la matriz P). _Estos arreglos luego son reducidos de manera concurrente para obtener los valores finales. 

_Nota: Al reducir el arreglo se probaron dos enfoques:_

* _Delegar la tarea de reducción al proceso “root” (idP = 0)_
* _Paralelizar la reducción del arreglo entre los T procesos._

_Calculando tiempos de ejecución se observó que cuando la matriz es relativamente chica y los procesos son pocos es menor el tiempo si la reducción la realiza el root. Sin embargo, la estrategia óptima en cuanto a escalabilidad y eficiencia es implementar una reducción de los elementos recopilados distribuyendo el trabajo entre los procesos._

##### OPENMP

Para la implementación en OpenMP se siguió una estrategia similar pero reescribiendo los algoritmos con las directivas específicas de la librería. Se utilizó _<span style="text-decoration:underline;">#pragma omp parallel</span>_ para distribuir el trabajo entre los T procesos. Luego para distribuir los bloques de matrices, lo que anteriormente llamamos _“distribuir el for externo”, _se utilizó la directiva _<span style="text-decoration:underline;">#pragma omp for</span>_. Estos bucles corresponden, cada uno, a una de las etapas descritas en el algoritmo secuencial. Finalmente se utilizó la directiva _<span style="text-decoration:underline;">reduction</span>_ dentro de los bucles principales para obtener la suma total de los elementos de la matriz P, para obtener el máximo de la matriz D y el mínimo de la matriz A.

##### MPI

Trabajar con MPI y en una arquitectura distribuida supuso un desafío mayor ya que no solo se debió repartir el trabajo entre los cores como en los programas anteriores, sino que además se debió garantizar que cada core tuviera acceso a los datos necesarios para realizar la operación.

Luego de un detenido análisis se observó que los procesos sólo necesitan conocer de las matrices A y D, los bloques que procesan, mientras que de las matrices B y C necesitan todos los elementos ya que los procesan todos. Además para los productos escalares, todos los procesos necesitan conocer MaxD, MinA, y PromP.

Para los elementos que son requeridos por todos los procesos, se delegó al proceso root (_idP = 0)_ realizar un _Broadcast _(comunicación colectiva).

Para las matrices que se reparten entre procesos (A y D) se realizó un _Scatter_, donde a cada proceso se le asigna una cantidad de elementos _equitativa_. En este punto surge la primera _limitación_ con respecto a los otros algoritmos: si la cantidad de procesos es tal que la cantidad de elementos no es perfectamente divisible por esta, el comportamiento es indeterminado. Sin embargo, esta limitación no es un problema para nuestro caso de uso específico, donde las matrices son de tamaño 1024, 2048 y 4096; con 4, 8 y 16 procesos.

Una vez dispersados los elementos de las matrices a cada proceso, se realizan las operaciones de cada etapa de manera independiente. Entre cada etapa es necesario realizar reducciones y comunicaciones únicamente para obtener las variables mínimo, máximo y promedio. El _Gather_ de los elementos dispersos sólo es necesario al final de la operación para obtener la matriz R y validar los resultados.

## Análisis de escalabilidad y conclusiones.

Para evaluar la estrategia de paralelización utilizada y la escalabilidad de la solución, se utilizarán las métricas vistas en clase (_speedup y eficiencia_). Para obtener los tiempos de ejecución necesarios para las métricas de cada programa, se utilizó el cluster provisto por la cátedra. A continuación se comparten los resultados obtenidos para los programas paralelos.

#### Tiempo de ejecución: Pthreads

<table>
  <tr>
   <td>
   </td>
   <td colspan="3" ><strong>Tamaño de problema (N)</strong>
   </td>
  </tr>
  <tr>
   <td><strong>Unidades de procesamiento</strong>
   </td>
   <td><strong>1024</strong>
   </td>
   <td><strong>2048</strong>
   </td>
   <td><strong>4096</strong>
   </td>
  </tr>
  <tr>
   <td><em>Secuencial</em>
   </td>
   <td><em>49.532s</em>
   </td>
   <td><em>410.882s</em>
   </td>
   <td>3289.593s
   </td>
  </tr>
  <tr>
   <td><strong>4</strong>
   </td>
   <td>12.298s
   </td>
   <td>105.354s
   </td>
   <td>833.612s
   </td>
  </tr>
  <tr>
   <td><strong>8</strong>
   </td>
   <td>6.167s
   </td>
   <td>52.149s
   </td>
   <td>413.710s
   </td>
  </tr>
</table>

#### Tiempo de ejecución: OpenMP

<table>
  <tr>
   <td>
   </td>
   <td colspan="3" ><strong>Tamaño de problema (N)</strong>
   </td>
  </tr>
  <tr>
   <td><strong>Unidades de procesamiento</strong>
   </td>
   <td><strong>1024</strong>
   </td>
   <td><strong>2048</strong>
   </td>
   <td><strong>4096</strong>
   </td>
  </tr>
  <tr>
   <td><em>Secuencial</em>
   </td>
   <td><em>49.532s</em>
   </td>
   <td><em>410.882s</em>
   </td>
   <td>3289.593s
   </td>
  </tr>
  <tr>
   <td><strong>4</strong>
   </td>
   <td>12.439s
   </td>
   <td>103.409s
   </td>
   <td>826.278s
   </td>
  </tr>
  <tr>
   <td><strong>8</strong>
   </td>
   <td>6.243s
   </td>
   <td>51.852s
   </td>
   <td>416.789s
   </td>
  </tr>
</table>

#### Tiempo de ejecución: MPI

<table>
  <tr>
   <td>
   </td>
   <td colspan="3" ><strong>Tamaño de problema (N)</strong>
   </td>
  </tr>
  <tr>
   <td><strong>Unidades de procesamiento</strong>
   </td>
   <td><strong>1024</strong>
   </td>
   <td><strong>2048</strong>
   </td>
   <td><strong>4096</strong>
   </td>
  </tr>
  <tr>
   <td><em>Secuencial</em>
   </td>
   <td><em>49.532s</em>
   </td>
   <td><em>410.882s</em>
   </td>
   <td>3289.593s
   </td>
  </tr>
  <tr>
   <td><strong>4</strong>
   </td>
   <td>12.762s
   </td>
   <td>104.427s
   </td>
   <td>834.802s
   </td>
  </tr>
  <tr>
   <td><strong>8</strong>
   </td>
   <td>6.551s
   </td>
   <td>52.039s
   </td>
   <td>411.802s
   </td>
  </tr>
  <tr>
   <td><strong>16</strong>
   </td>
   <td>3.499s
   </td>
   <td>26.960s
   </td>
   <td>213.770s
   </td>
  </tr>
</table>

### SPEEDUP 

Es una métrica de rendimiento que muestra el incremento relativo obtenido al ejecutar un programa luego de aplicar mejoras. A nuestro programa inicial secuencial se le realizaron mejoras significativas por medio de la paralelización.

#### SPEEDUP para Pthreads

<table>
  <tr>
   <td>
   </td>
   <td colspan="3" ><strong>Tamaño de problema (N)</strong>
   </td>
  </tr>
  <tr>
   <td><strong>Unidades de procesamiento</strong>
   </td>
   <td><strong>1024</strong>
   </td>
   <td><strong>2048</strong>
   </td>
   <td><strong>4096</strong>
   </td>
  </tr>
  <tr>
   <td><strong>4</strong>
   </td>
   <td>3.88
   </td>
   <td>3.9
   </td>
   <td>3.95
   </td>
  </tr>
  <tr>
   <td><strong>8</strong>
   </td>
   <td>7.93
   </td>
   <td>7.88
   </td>
   <td>7.95
   </td>
  </tr>
</table>

#### SPEEDUP para OpenMP

<table>
  <tr>
   <td>
   </td>
   <td colspan="3" ><strong>Tamaño de problema (N)</strong>
   </td>
  </tr>
  <tr>
   <td><strong>Unidades de procesamiento</strong>
   </td>
   <td><strong>1024</strong>
   </td>
   <td><strong>2048</strong>
   </td>
   <td><strong>4096</strong>
   </td>
  </tr>
  <tr>
   <td><strong>4</strong>
   </td>
   <td>3.98
   </td>
   <td>3.97
   </td>
   <td>3.98
   </td>
  </tr>
  <tr>
   <td><strong>8</strong>
   </td>
   <td>7.93
   </td>
   <td>7.92
   </td>
   <td>7.89
   </td>
  </tr>
</table>

#### SPEEDUP para MPI

<table>
  <tr>
   <td>
   </td>
   <td colspan="3" ><strong>Tamaño de problema (N)</strong>
   </td>
  </tr>
  <tr>
   <td><strong>Unidades de procesamiento</strong>
   </td>
   <td><strong>1024</strong>
   </td>
   <td><strong>2048</strong>
   </td>
   <td><strong>4096</strong>
   </td>
  </tr>
  <tr>
   <td><strong>4</strong>
   </td>
   <td>3.88
   </td>
   <td>3.93
   </td>
   <td>3.94
   </td>
  </tr>
  <tr>
   <td><strong>8</strong>
   </td>
   <td>7.56
   </td>
   <td>7.90
   </td>
   <td>7.99
   </td>
  </tr>
  <tr>
   <td><strong>16</strong>
   </td>
   <td>14.15
   </td>
   <td>15.24
   </td>
   <td>15.43
   </td>
  </tr>
</table>

### EFICIENCIA

Medida de rendimiento que indica el porcentaje de tiempo en el que las unidades de procesamiento realizan trabajo útil.

#### EFICIENCIA para Pthreads

<table>
  <tr>
   <td>
   </td>
   <td colspan="3" ><strong>Tamaño de problema (N)</strong>
   </td>
  </tr>
  <tr>
   <td><strong>Unidades de procesamiento</strong>
   </td>
   <td><strong>1024</strong>
   </td>
   <td><strong>2048</strong>
   </td>
   <td><strong>4096</strong>
   </td>
  </tr>
  <tr>
   <td><strong>4</strong>
   </td>
   <td>0.97
   </td>
   <td>0.975
   </td>
   <td>0.988
   </td>
  </tr>
  <tr>
   <td><strong>8</strong>
   </td>
   <td>0.998
   </td>
   <td>0.985
   </td>
   <td>0.994
   </td>
  </tr>
</table>

#### EFICIENCIA para OpenMP

<table>
  <tr>
   <td>
   </td>
   <td colspan="3" ><strong>Tamaño de problema (N)</strong>
   </td>
  </tr>
  <tr>
   <td><strong>Unidades de procesamiento</strong>
   </td>
   <td><strong>1024</strong>
   </td>
   <td><strong>2048</strong>
   </td>
   <td><strong>4096</strong>
   </td>
  </tr>
  <tr>
   <td><strong>4</strong>
   </td>
   <td>0.995
   </td>
   <td>0.993
   </td>
   <td>0.995
   </td>
  </tr>
  <tr>
   <td><strong>8</strong>
   </td>
   <td>0.991
   </td>
   <td>0.99
   </td>
   <td>0.99
   </td>
  </tr>
</table>

#### EFICIENCIA para MPI

<table>
  <tr>
   <td>
   </td>
   <td colspan="3" ><strong>Tamaño de problema (N)</strong>
   </td>
  </tr>
  <tr>
   <td><strong>Unidades de procesamiento</strong>
   </td>
   <td><strong>1024</strong>
   </td>
   <td><strong>2048</strong>
   </td>
   <td><strong>4096</strong>
   </td>
  </tr>
  <tr>
   <td><strong>4</strong>
   </td>
   <td>0.97
   </td>
   <td>0.98
   </td>
   <td>0.98
   </td>
  </tr>
  <tr>
   <td><strong>8</strong>
   </td>
   <td>0.94
   </td>
   <td>0.987
   </td>
   <td>0.999
   </td>
  </tr>
  <tr>
   <td><strong>16</strong>
   </td>
   <td>0.88
   </td>
   <td>0.952
   </td>
   <td>0.962
   </td>
  </tr>
</table>

### ESCALABILIDAD

Un programa paralelo es escalable si mantiene constante la eficiencia al aumentar el número de unidades de procesamiento aumentando también el tamaño del problema.

#### Observando los resultados obtenidos para Pthreads y OpenMP

El programa resulta fuertemente escalable ya que mantiene el nivel de eficiencia relativamente constante a medida que el tamaño del problema aumenta.

#### Observando los resultados obtenidos para MPI

Se observa para el caso de 8 y 16 unidades de procesamiento una disminución de la eficiencia cuando el tamaño de la matriz es de 1024 y 2048, esto se debe posiblemente a que el tamaño del bloque de la matriz se ve reducido, lo que afecta la localidad de las operaciones matriciales. 

Además, para 16 unidades de procedimiento la eficiencia se ve reducida con respecto a los otros casos, un factor que podría perjudicar a la escalabilidad es la necesidad de aumentar las comunicaciones entre las unidades de procesamiento.
