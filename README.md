# filtrado_lecturas_cooledKLL

Este proyecto consiste en la implementación del sketch Cooled-KLL y su aplicación para filtrar aquellos k-mers que se encuentren en cuantiles extremos de la distribución y que puedan ser fruto de errores de secuenciación, los cuales pueden sesgar y disminuir la precisión de los resultados en los estudios de secuenciación masiva. También se realizan múltiples experimentos para evaluar el rendimiento del sketch en base al dataset de 50 genomas de Escherichia Coli.

# Filtrado de lecturas

Para ejecutar el filtrado de lecturas ejecute los siguientes comandos:

```bash
g++ -std=c++20 -o filtrar_kmers source/filtrar_kmers.cpp
./filtrar_kmers <folder_file> <save_file> <k-mers_length> <N_buckets> <B_capacity> <C_size> <l_quantile> <u_quantile>
```
Donde:

**<folder_file>:** ruta a la carpeta con los archivos genomicos de tipo FASTA.
**<save_file>:** ruta al archivo donde se guardan las estadisticas del filtro de datos.
**<k-mers_length>:** largo de los k-mers.
**<N_buckets>:** número de buckets en el hot filter del sketch.
**<B_capacity>:** número de entradas que tiene cada bucket.
**<C_size>:** número de elementos en el compactador más grande en el KLL clasico.
**<l_quantile>:** cuantil inferior para filtrar los datos.
**<u_quantile>:** cuantil superior para filtrar los datos.

<small>**la carpeta indicada por <folder_file> debe contener una serie de archivos de tipo FASTA con datos genomicos.**</small>

# Estimación de la distribución de los datos

A continuación se detallan los pasos a realizar para estimar la distribución de abundancia de los k-mers obtenidos a partir de un conjunto de lecturas y obtener un CSV con datos sobre la distribución estimada y la real.

## 1. Sin archivos previos:

Compilar el archivo **LeerYEstimarDistribucion.cpp** y ejecutar con los siguientes comandos:
```bash
g++ -std=c++20 -o estimar_distribucion source/estimar_distribucion.cpp
./estimar_distribucion <kmers_file> <k-mers_length> <distribution> <N_buckets> <B_capacity> <C_size>
```
Donde:

**<folder_file>:** ruta a la carpeta con los datos genomicos.
**<k-mers_length>:** largo de los k-mers.
**\<distribution>:** define en base a que variable se calcula la distribución: 0 = distribucion de k-mers | 1 = distribución de frecuencias.
**<N_buckets>:** número de bloques en el hot filter del sketch.
**<B_capacity>:** número de entradas de cada bloque.
**<C_size>:** número de elementos en el compactor más grande en la parte KLL clasico.

<small>**la carpeta indicada por <folder_file> debe contener una serie de archivos de tipo FASTA con datos genomicos.**</small>

Al ejecutar estos comandos el programa leera todos los archivos FASTA en el directorio indicado, obtendra los k-mers a partir de las lecturas leidas, estimara la distribución a partir de un Cooled-KLL construido con los parametros ingresados y obtendra la distribución real mediante un método determinista. Finalmente, en la carpeta **data/frequency_distribution** se creara una carpeta, en caso de no existir, con un nombre basado en los parámetros de construcción del sketch y allí se guardara una archivo CSV con los datos de las distribuciones obtenidas y otro con la memoria ocupada por los distintos métodos.

## 2. Con los k-mers ya almacenados:

Para almacenar los k-mers en un archivo y poder ejecutar los experimentos con distintas configuraciones de sketch más rápidamente puedes seguir los siguientes pasos:

1. Compilar y ejecutar **leer_kmers.cpp** para obtener los k-mers con sus frecuencias y almacenarlos en un CSV.

    ```bash
    g++ -std=c++20 leer_kmers source/leer_kmers.cpp
    ./leer_kmers <folder_url> <k>
    ```
    **<folder_url>:** path to the folder where FASTA files are located.
    **<k>:** length of the kmer.

    Luego de la ejecución, en la carpeta **data/kmers** se creara un archivo CSV con los k-mers y sus frecuencias presentes en las lecturas leídas y sus frecuencias.

2. Compilar y ejecutar **estimar_distribucion.cpp** para obtener los resultados de las distribuciones.

    ```bash
    g++ -std=c++20 estimar_distribucion source/estimar_distribucion.cpp
    ./estimar_distribucion <kmers_file> <k-mers_length> <distribution> <N_buckets> <B_capacity> <C_size>
    ```

    Donde:

    **<kmers_file>:** ruta hacia el archivo con los k-mers.
    **<k-mers_length>:** largo de los k-mers en el archivo.
    **\<distribution>:** define en base a que variable se calcula la distribución: 0 = distribucion de k-mers | 1 = distribución de frecuencias.
    **<N_buckets>:** número de bloques en el hot filter del sketch.
    **<B_capacity>:** número de entradas de cada bloque.
    **<C_size>:** número de elementos en el compactor más grande en la parte KLL clasico.

    El resultado es el mismo que en la ejecución sin archivos previos.

En ambos métodos loas archivos que contienen la información de los resultados de la estimación de la distribución serán almacenados en una carpeta con nombre basado en los parámetros de construcción del sketch, y el nombre del archivo comenzara con el tamaño de k-mer y seguira con **mers_distribution** o **mers_memory** dependiendo de lo que almacene. Esta división se realizo para obtener distintas graficas según la configuración del sketch.

Adicionalmente, si se crean un conjunto de archivos almacenando los k-mers y sus frecuencias en la carpeta **data/kmers** mediante el archivo **leer_kmers.cpp**, es posible ejecutar la estimación de las distribuciones para los distintos tamaños de k-mer usando el archivo bash **estimar_distribuciones.sh**. Para eso es necesario compilar el archivo **estimar_distribucion.cpp** dentro de una carpeta **bin/**, para eso primero cree la carpeta, compile el archivo y luego ejecute el siguiente comando:

```bash
./estimar_distribuciones.sh <N_buckets> <B_capacity> <C_size> <distribution>
```

Donde:

**<N_buckets>:** número de bloques en el hot filter del sketch.
**<B_capacity>:** número de entradas de cada bloque.
**<C_size>:** número de elementos en el compactor más grande en la parte KLL clasico.
**\<distribution>:** define en base a que variable se calcula la distribución: 0 = distribucion de k-mers | 1 = distribución de frecuencias.

# Experimentación

Se realizaron dos experimentos a las dos estructuras deterministas y a la estructura estocastica, uno sobre el tamaño de construcción de las estructuras y otro sobre el tiempo usado para realizar las consultas quantile($\delta$) y rank($x$).

## Experimento de tiempo de construcción

Este experimento mide el tiempo que tardan en construir cada una de las soluciones planteadas.

```bash
g++ -std=c++20 -o uhr_construccion source/uhr_construccion.cpp
./uhr_construccion <filename> <RUNS> <METHOD>
```
Donde:

**\<filename>:** ruta y nombre del archivo donde los resultados del experimento serán escritos (con extension .csv).
**\<RUNS>:** número de ejecuciones por caso de prueba: debería ser >= 32.
**\<METHOD>:** 1 = vector plano | 2 = vector comprimido | 3 = sketch.

Si se quiere modificar la configuración del sketch, es necesario modificarla manualmente en la linea 168 y 169 del código.

## Experimento de tiempo de consultas

Este experimento mide cuanto demoran en realizar las consultas quantile($\delta$) y rank(x) cada una de las soluciones propuestas.

```bash
g++ -std=c++20 -o uhr_quantile_rank uhr_quantile_rank.cpp
./uhr_quantile_rank <filename> <RUNS> <METHOD> <k>
```

Donde:

**\<filename>:** ruta y nombre del archivo donde los resultados del experimento serán escritos (con extension .csv).
**\<RUNS>:** número de ejecuciones por caso de prueba: debería ser >= 32.
**\<METHOD>:** 1 = vector plano | 2 = vector comprimido | 3 = sketch.
**\<k>**: largo del k-mer a utilizar.

Si se quiere modificar la configuración del sketch, es necesario modificarla manualmente en la linea 132 y 133 del código.

Adicionalmente, para automatizar los experimentos se puede usar el archivo  **run_experiments.sh**, el cual ejecuta los experimentos de construcción y consulta para las tres soluciones distintas y en caso de los experimentos de consulta, los ejecuta para k $\in [3, 6, 9, ..., 27, 30]$. Para usar este archivo es necesario tener la carpeta **data/experiments/** creada.

# Obtención de graficos

Una vez que se obtienen los archivos CSV con los resultados de la distribución de frecuencias estimada y real, podemos continuar con la obtención de graficos.

## Instalar las librerias necesarias

```bash
pip install pandas matplotlib numpy adjustText
```

## Obtención de graficos
### 1. Gráficos distribución real vs estimada

Para obtener los graficos de distribucion real vs estimada, es necesario que modifique la variable **NOMBRE_DIRECTORIO** con el nombre del directorio en el que se encuentran los CSV con los resultados de la estimación de las distribuciones. Luego, ejecute el archivo **graficos_distribucion.py** con el siguiente comando:

```bash
python source/graficos_distribucion.py
```

En la carpeta **data/images/** se generara un conjunto de imagenes con las distribuciones para cada uno de los valores de k usados.

### 2. Gráficos error relativo y absoluto en la estimación de las distribuciones

Para obtener los gráficos de error para cada una de las distribuciones obtenidas con una configuración especifica del sketch, modifique la variable **SKETCH_SETTINGS** en el archivo **graficos_error_kmer.py** con el nombre de la carpeta encontrada en **data/frequency_distribution/** de la que se quieren obtener los gráficos. Luego ejecute:

```bash
python source/graficos_error_kmer.py
```

A continuación, se generaran los graficos en la carpeta **data/images/**.

### 3. Gráficos tiempos de consulta

Luego de haber realizada la experimentación de tiempo de consulta podemos obtener un gráfico comparativo de las distintas estructuras, para esto es necesario haber realizado la experimentación con las tres estructuras para distintos largos de k-mers. Para obtener los gráficos ejecute el siguiente comando:

```bash
python source/graficos_consultas.py
```

Luego de ejecutado el comando se podran encontrar los gráficos en la carpeta **data/images/graficos_consultas/**.

### 4. Gráfico tiempo de construcción

Luego de realizada la experimentación de tiempo de construcción, podemos obtener un gráfico comparativo de las distintas estructuras, para ello es necesario haber ejecutado la experimentación con las tres estructuras. Para
obtener los gráficos ejecute el siguiente comando:

```bash
python source/graficos_exp_construccion.py
```

El gráfico resultante se encontrara en la carpeta **data/images/**.

### 5. Tabla estadísticas de filtrado

Para tabular las estadísticas del filtrado de lecturas es necesario modificar en el archivo **tabla_filtrado_lecturas.py** la variable **ARCHIVO_ENTRADA** con la ruta hacia el CSV con las estadísticas del filtrado obtenidas previamente y **ARCHIVO_SALIDA** con la ruta y nombre del archivo donde se quiere almacenar la tabla. Luego, ejecutar el siguiente comando:

```bash
python source/tabla_filtrado_lecturas.py
```

### 6. Tabla uso de memoria

Para tabular los datos de la memoria usada por cada una de las soluciones es necesario modificar en el archivo **tabla_uso_memoria.py** la variable **SKETCH_SETTINGS** con el nombre de la carpeta que contiene los resultados de la estimación de distribución para una configuración de sketch especifica. Para obtener la tabla, ejecute el siguiente comando:

```bash
python source/tabla_uso_memoria.py
```

La tabla generada se podrá encontrar en la carpeta **data/images/SKETCH_SETTINGS**



    