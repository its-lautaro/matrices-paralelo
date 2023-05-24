# Sistemas Distribuidos y Paralelos
## Trabajo de promoción

### Calculo del tamaño óptimo del bloque de memoria

| N / BS   | 32      | 64      | 128         | 256         | 512        | 1024    | 2048    | 4096 |
| -------- | ------- | ------- | ----------- | ----------- | ---------- | ------- | ------- | ---- |
| **1024** | 13.207  | 13.038  | 12.873      | 12.797      | **12.761** | -       | -       | -    |
| **2048** | 105.212 | 104.202 | 102.841     | **102.329** | 105.564    | 104.719 | -       | -    |
| **4096** | 839.167 | 833.465 | **823.959** | 855.641     | 846.233    | 838.107 | 833.543 | -    |


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

### Algoritmo secuencial
Para desarrollar el algoritmo secuencial seguimos dos alternativas. En una primera evaluamos el tiempo realizando cada operación en funciones separadas, sin considerar la dependencia de datos.
| 1024   | 2048    | 4096     |
| ------ | ------- | -------- |
| 51.278 | 409.614 | 3294.489 |

Para la segunda versión del algoritmo secuencial, intentamos aprovechar la dependencia de datos realizando en una misma iteración las operaciones que utilizaban los mismos datos:
- El producto A.B.C y el producto D.C.B, que utilizan las matrices C y B
- Junto con el producto de las matrices, también calculamos el minimo y el maximo en las matrices A y D respectivamente.
- Luego al calcular los elementos de la matriz P también realizamos la suma de estos elementos para el cálculo del promedio
- Finalmente realizamos el cálculo de R usando el promedio y la matriz P calculados anteriormente. Esta operación solo puede realizarse una vez calculado el promedio por lo que debemos volver a procesar la matriz P.
  
| 1024   | 2048    | 4096     |
| ------ | ------- | -------- |
| 50.226 | 408.717 | 3289.593 |

