# Sistemas Distribuidos y Paralelos
## Trabajo de promoción

### Calculo del tamaño óptimo del bloque de memoria

`sbatch script_sec.sh 1024 128 & sbatch script_sec.sh 1024 256 & sbatch script_sec.sh 1024 512`

|      | 128    | 256    | 512    | 1024   | 2048 |
| ---- | ------ | ------ | ------ | ------ | ---- |
| 1024 | 12.899 | 12.800 | 12.765 |  | -    |
| 2048 | 80.926 | 80.248 | 80.985 |        |      |
| 4096 |        |        |        |        |      |
