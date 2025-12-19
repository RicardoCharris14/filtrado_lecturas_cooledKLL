#!/bin/bash

# --- Configuración de parámetros ---

# El script espera 3 argumentos: los últimos 3 parámetros para el ejecutable.
if [ "$#" -ne 4 ]; then
    echo "Correct usage: $0 <N_buckets> <B_capacity> <C_size> <distribution>"
    echo "<N_buckets>: number of buckets in the hot filter part of the sketch."
    echo "<B_capacity>: number of entries a bucket have."
    echo "<C_size>: number of elements of the largest compactor in the classic kll part."
    echo "<distribution>: 1 = abundance distribution | 0 = kmers distribution."
    exit 1
fi

# Asignar los argumentos de entrada a variables
PARAM3=$1  # e.g., 10
PARAM4=$2  # e.g., 5
PARAM5=$3  # e.g., 20

# Ubicación del ejecutable (usando el formato de path para Windows en entornos WSL/Git Bash)
EXECUTABLE="./bin/estimar_distribucion.exe"
# Directorio donde se encuentran los archivos CSV
DATA_DIR="./data/kmers"

# Parámetro fijo que va después del tamaño del k-mer
FIXED_ARG=$4

echo "--- Iniciando experimentos ---"
echo "Ejecutable: $EXECUTABLE"
echo "Directorio de datos: $DATA_DIR"
echo "Parámetros fijos: $FIXED_ARG $PARAM3 $PARAM4 $PARAM5"
echo "------------------------------"

# Iterar sobre todos los archivos que coinciden con el patrón *mers_frequency.csv
# El patrón se ajusta a lo que se ve en la imagen: 5mers_frequency.csv, 21mers_frequency.csv, etc.
find "$DATA_DIR" -name "*mers_frequency.csv" -print0 | while IFS= read -r -d $'\0' csv_file; do

    # 1. Obtener el nombre base del archivo (ej: 5mers_frequency.csv)
    filename=$(basename "$csv_file")

    # 2. Extraer el número inicial (el tamaño del k-mer)
    # Usa 'expr' con RegEx para buscar el patrón 'NUMERO'mers
    # El resultado se asigna a kmer_size.
    # Ejemplo: '5mers_frequency.csv' -> '5'
    kmer_size=$(expr "$filename" : '\([0-9]\+\)mers')

    # Verificar que se haya extraído un tamaño de k-mer válido
    if [ -z "$kmer_size" ]; then
        echo "Advertencia: No se pudo extraer el tamaño del k-mer de: $filename. Saltando..."
        continue
    fi

    # 3. Construir y mostrar el comando a ejecutar
    COMMAND="$EXECUTABLE $csv_file $kmer_size $FIXED_ARG $PARAM3 $PARAM4 $PARAM5 < /dev/null"
    echo "--> Ejecutando: $COMMAND"

    # 4. Ejecutar el comando
    # El uso de 'cmd.exe /C' es a menudo necesario si estás ejecutando un .exe de Windows
    # desde un entorno Linux-like (como WSL o Git Bash)
    "$EXECUTABLE" "$csv_file" "$kmer_size" "$FIXED_ARG" "$PARAM3" "$PARAM4" "$PARAM5" < /dev/null
    
    # 5. Comprobar si la ejecución fue exitosa
    if [ $? -ne 0 ]; then
        echo "!!! Error en la ejecución para: $filename (código de salida: $?)."
    else
        echo "--- Éxito para: $filename"
    fi

    echo "" # Línea vacía para separar las salidas

done

echo "--- Todos los experimentos han finalizado ---"