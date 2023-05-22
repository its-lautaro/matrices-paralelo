# Sistemas Distribuidos y Paralelos
## Trabajo de promoción

### Calculo del tamaño óptimo del bloque de memoria

`sbatch script_sec.sh 1024 128 & sbatch script_sec.sh 1024 256 & sbatch script_sec.sh 1024 512`

>Promedio de 3 ejecuciones

|      | 128     | 256     | 512     | 1024    | 2048    |
| ---- | ------- | ------- | ------- | ------- | ------- |
| 1024 | 12.886  | 12.792  | 12.751  | 12.935  | -       |
| 2048 | 102.952 | 102.215 | 105.622 | 104.657 | 103.543 |
| 4096 |         |         |         |         |         |
