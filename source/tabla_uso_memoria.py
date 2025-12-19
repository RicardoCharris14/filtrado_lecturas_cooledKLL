import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import re

# --- CONFIGURACIÓN ---
# Ajusta la ruta a tu carpeta
SKETCH_SETTINGS = 'NB_100_BC_10_CS_100'

CARPETA_DATOS = 'data/frequency_distribution/'+SKETCH_SETTINGS
CARPETA_SALIDA = 'data/images/'+SKETCH_SETTINGS

# Columnas originales del CSV
COLUMNAS_DESEADAS = [
    'elements', 
    'unique_elements', 
    'sketch_memory', 
    'vector_memory', 
    'compressed_vector_memory'
]

# Columnas que representan memoria (para dividirlas por 1024)
COLUMNAS_MEMORIA = [
    'sketch_memory', 
    'vector_memory', 
    'compressed_vector_memory'
]

# Nombres bonitos para la tabla (Agregamos "(KB)")
NOMBRES_COLUMNAS = [
    'K-mer', 
    'Elements', 
    'Unique Elem.', 
    'Sketch Mem (KB)', 
    'Vector Mem (KB)', 
    'Comp. Vector Mem (KB)'
]

def generar_tabla_memoria_kb():
    print(f"--- Buscando archivos de memoria en: {CARPETA_DATOS} ---")
    
    if not os.path.exists(CARPETA_DATOS):
        print("[ERROR] La carpeta de datos no existe.")
        return

    patron = os.path.join(CARPETA_DATOS, '*mers_memory.csv')
    archivos = glob.glob(patron)
    
    if not archivos:
        print("[ERROR] No se encontraron archivos '*mers_memory.csv'.")
        return

    lista_filas = []

    for ruta in archivos:
        nombre_archivo = os.path.basename(ruta)
        
        # 1. Extraer K
        match = re.search(r'\d+', nombre_archivo)
        if match:
            k_val = int(match.group())
        else:
            continue
            
        try:
            # 2. Leer CSV
            df = pd.read_csv(ruta)
            if df.empty: continue
            
            # 3. Seleccionar columnas
            cols_existentes = [c for c in COLUMNAS_DESEADAS if c in df.columns]
            fila = df[cols_existentes].iloc[0].to_dict()
            
            # 4. Agregar K
            fila['K'] = k_val
            
            lista_filas.append(fila)
            
        except Exception as e:
            print(f"[ERROR] Leyendo {nombre_archivo}: {e}")

    if not lista_filas:
        print("No se extrajeron datos válidos.")
        return

    # --- PROCESAMIENTO ---
    df_final = pd.DataFrame(lista_filas)
    
    # 1. Ordenar por K
    df_final = df_final.sort_values(by='K')
    
    # 2. CONVERSIÓN A KB (Dividir por 1024)
    for col in COLUMNAS_MEMORIA:
        if col in df_final.columns:
            df_final[col] = df_final[col] / 1024
            # Redondear a 2 decimales para que no salga 12.3421515
            df_final[col] = df_final[col].round(2)

    # 3. Reordenar columnas (K primero)
    cols = ['K'] + [c for c in COLUMNAS_DESEADAS if c in df_final.columns]
    df_final = df_final[cols]

    # --- GRAFICAR TABLA ---
    alto_fig = max(3, len(df_final) * 0.5 + 1)
    fig, ax = plt.subplots(figsize=(13, alto_fig)) # Un poco más ancha para que quepa "(KB)"
    
    ax.axis('off')
    ax.axis('tight')
    
    tabla = ax.table(
        cellText=df_final.values,
        colLabels=NOMBRES_COLUMNAS,
        cellLoc='center',
        loc='center'
    )
    
    tabla.auto_set_font_size(False)
    tabla.set_fontsize(11)
    tabla.scale(1.2, 2)
    
    # Estilo de encabezado y filas
    for (row, col), cell in tabla.get_celld().items():
        if row == 0:
            cell.set_text_props(weight='bold', color='white')
            cell.set_facecolor('#40466e')
        elif row % 2 == 0:
            cell.set_facecolor('#f1f1f2')

    plt.title(f"Resumen de Memoria por K-mer (KB)\n({os.path.basename(CARPETA_DATOS)})", 
              weight='bold', fontsize=14, pad=20)

    if not os.path.exists(CARPETA_SALIDA):
        os.makedirs(CARPETA_SALIDA)
        
    ruta_salida = os.path.join(CARPETA_SALIDA, 'tabla_resumen_memoria.png')
    plt.savefig(ruta_salida, bbox_inches='tight', dpi=150)
    plt.close(fig)
    
    print(f"Tabla (en KB) guardada en: {ruta_salida}")
    print("\n--- Datos Procesados (KB) ---")
    print(df_final.to_string(index=False))

if __name__ == "__main__":
    generar_tabla_memoria_kb()