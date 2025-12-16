#!/bin/bash

# --- Configuración ---
EXE_CONSTRUCCION="./bin/uhr_construccion.exe"
EXE_CONSULTAS="./bin/uhr_quantile_rank.exe"

OUTPUT_DIR="data/experiments"
BASE_NAME_CONST="construccion"
BASE_NAME_QUERY="consultas"

# Número de corridas (Ajustar según necesidad, ej: 12 o 32)
RUNS=12 

mkdir -p "$OUTPUT_DIR"

echo "==========================================="
echo "   INICIANDO BATERÍA DE EXPERIMENTOS"
echo "==========================================="

# Bucle principal: Métodos (1, 2, 3)
for METHOD in 1 2 3; do

    # Definir sufijo base según el método
    case $METHOD in
        1) SUFFIX_BASE="plain_v" ;;
        2) SUFFIX_BASE="compressed_v" ;;
        3) SUFFIX_BASE="sketch" ;;
    esac

    echo ""
    echo ">>> PROCESANDO METODO $METHOD ($SUFFIX_BASE)"
    echo "-------------------------------------------"

    # -------------------------------------------------------
    # 1. EXPERIMENTO DE CONSTRUCCIÓN (Se ejecuta una vez por Método)
    # -------------------------------------------------------
    # Nota: Asumimos que la construcción NO depende de k, por eso va fuera del loop de k.
    
    FILE_CONST="${OUTPUT_DIR}/${BASE_NAME_CONST}_${SUFFIX_BASE}.csv"
    
    if [ -f "$FILE_CONST" ]; then
        echo " [Construccion] Aviso: $FILE_CONST ya existe (se podría sobrescribir)."
    fi

    echo " -> Ejecutando Construcción..."
    "$EXE_CONSTRUCCION" "$FILE_CONST" "$RUNS" "$METHOD" < /dev/null

    if [ $? -ne 0 ]; then
        echo "    [ERROR] Falló la construcción. Saltando consultas para este método."
        continue 
    else
        echo "    [OK] Construcción terminada."
    fi

    # -------------------------------------------------------
    # 2. EXPERIMENTO DE CONSULTAS (Iterando k de 3 a 30, salto de 3)
    # -------------------------------------------------------
    echo " -> Iniciando barrido de k (3..30)..."

    # seq 3 3 30 genera: 3, 6, 9, 12, ..., 30
    for k in $(seq 3 3 30); do
        
        # Generamos un nombre único por k para no perder datos
        # Ej: data/experiments/consultas_plain_v_k3.csv
        FILE_QUERY="${OUTPUT_DIR}/${BASE_NAME_QUERY}_${SUFFIX_BASE}_k${k}.csv"

        echo "    Running k=$k -> $FILE_QUERY"

        # Ejecutamos pasando k como 4to parámetro
        # Estructura: <EXE> <FILE> <RUNS> <METHOD> <K>
        "$EXE_CONSULTAS" "$FILE_QUERY" "$RUNS" "$METHOD" "$k" < /dev/null

        if [ $? -ne 0 ]; then
            echo "    [ERROR] Falló k=$k"
        fi

    done

done

echo ""
echo "==========================================="
echo "   TODOS LOS EXPERIMENTOS FINALIZADOS"
echo "==========================================="