import pandas as pd
import matplotlib.pyplot as plt
import os

# --- CONFIGURACIÓN ---
# Ajusta la ruta si es necesario. Basado en tu imagen parece ser:
CARPETA_DATOS = 'data/experiments'
CARPETA_SALIDA = 'data/images'
# Diccionario con los archivos y el nombre que quieres que aparezca en la leyenda
ARCHIVOS = {
    'construccion_plain_v.csv': 'Plain Vector',
    'construccion_compressed_v.csv': 'Compressed Vector',
    'construccion_sketch.csv': 'Sketch'
}

if not os.path.exists(CARPETA_SALIDA):
    os.makedirs(CARPETA_SALIDA, exist_ok=True)
    print(f"Carpeta creada: {CARPETA_SALIDA}")

def graficar_tiempo_construccion():
    print(f"--- Iniciando gráfico de Tiempos de Construcción ---")
    
    # Crear la figura
    fig, ax = plt.subplots(figsize=(10, 6))

    archivos_encontrados = False

    for nombre_archivo, etiqueta in ARCHIVOS.items():
        ruta_completa = os.path.join(CARPETA_DATOS, nombre_archivo)
        
        if not os.path.exists(ruta_completa):
            print(f"[AVISO] No se encontró: {nombre_archivo}")
            continue

        try:
            # Leer CSV
            # skipinitialspace=True ayuda si hay espacios después de las comas
            df = pd.read_csv(ruta_completa)

            df['t_mean'] = df['t_mean']/1e9
            df['t_stdev'] = df['t_stdev']/1e9

            # --- GRAFICAR CON BARRAS DE ERROR ---
            # x = n (largo del kmer)
            # y = t_mean (tiempo promedio)
            # yerr = t_stdev (desviación estándar)
            # capsize = tamaño de las "tapitas" en las barras de error
            ax.errorbar(
                df['n'], 
                df['t_mean'], 
                yerr=df['t_stdev'], 
                label=etiqueta, 
                marker='o',       # Puntos en cada dato
                capsize=5,        # Ancho de los topes de la barra de error
                linewidth=2,      # Grosor de la línea
                alpha=0.8         # Transparencia leve
            )
            
            print(f"Procesado: {nombre_archivo}")
            archivos_encontrados = True

        except Exception as e:
            print(f"[ERROR] Falló al leer {nombre_archivo}: {e}")

    if not archivos_encontrados:
        print("No se pudieron cargar datos. Verifica la ruta CARPETA_DATOS.")
        return

    # --- FORMATO DEL GRÁFICO ---
    ax.set_title("Largo del k-mer vs Tiempo de construcción", fontsize=14, fontweight='bold')
    ax.set_xlabel("K (Largo del K-mer)", fontsize=12)
    ax.set_ylabel("Tiempo Promedio (Segundos)", fontsize=12)
    
    # Grid para facilitar la lectura
    ax.grid(True, linestyle='--', alpha=0.6)
    
    # Leyenda
    ax.legend(fontsize=10)

    # Opcional: Si los tiempos varían demasiado (ej: de 100 a 1.000.000), 
    # descomenta la siguiente línea para usar escala logarítmica:
    # ax.set_yscale('log')

    plt.tight_layout()
    
    # Guardar y mostrar
    ruta_salida = os.path.join(CARPETA_SALIDA, 'grafico_tiempos_construccion.png')
    plt.savefig(ruta_salida, dpi=150)
    print(f"Gráfico guardado en: {ruta_salida}")
    plt.show()

if __name__ == "__main__":
    graficar_tiempo_construccion()