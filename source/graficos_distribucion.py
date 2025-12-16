import pandas as pd
import matplotlib.pyplot as plt
import os
import glob

# --- CONFIGURACIÓN ---
RUTA_BASE = 'data/frequency_distribution' 

# El directorio específico que quieres procesar
NOMBRE_DIRECTORIO = 'NB_20_BC_5_CS_20'

# Carpeta de salida
CARPETA_SALIDA = 'data/images/'+NOMBRE_DIRECTORIO

# Opción extra: Escala logarítmica
# Cambia a True si tienes frecuencias muy altas (ej: 10^6) y bajas a la vez
USAR_LOG_SCALE = False 

def procesar_grafico_unico():
    ruta_objetivo = os.path.join(RUTA_BASE, NOMBRE_DIRECTORIO)
    
    print(f"--- Iniciando análisis comparativo ---")
    print(f"Directorio: {NOMBRE_DIRECTORIO}")

    if not os.path.exists(ruta_objetivo):
        print(f"[ERROR] No existe el directorio: {ruta_objetivo}")
        return

    if not os.path.exists(CARPETA_SALIDA):
        os.makedirs(CARPETA_SALIDA)

    patron = os.path.join(ruta_objetivo, '*mers_distribution.csv')
    archivos = glob.glob(patron)
    archivos.sort()

    if not archivos:
        print("[AVISO] No se encontraron archivos CSV.")
        return

    count = 0
    for archivo_path in archivos:
        try:
            nombre_archivo = os.path.basename(archivo_path)
            label = nombre_archivo.replace('.csv', '') 
            
            df = pd.read_csv(archivo_path)
            if df.empty: continue

            # --- GENERAR GRÁFICO (1 SOLO AXIS) ---
            # Creamos una sola figura con un solo eje (ax)
            fig, ax = plt.subplots(figsize=(10, 6))
            
            titulo = f"Comparación Frecuencias: {label}\n({NOMBRE_DIRECTORIO})"
            ax.set_title(titulo, fontsize=12, fontweight='bold')

            # 1. Línea de Frecuencia REAL (Azul sólida)
            ax.plot(df['quantile'], df['real_quantile'], 
                    color='tab:blue', label='Real', linewidth=2, alpha=0.8)

            # 2. Línea de Frecuencia ESTIMADA (Naranja punteada)
            # Usamos linestyle='--' para que se vea la azul debajo si coinciden
            ax.plot(df['quantile'], df['estimated_quantile'], 
                    color='tab:orange', label='Estimada', linewidth=2, linestyle='--')

            ax.set_xlabel("Quantile")
            ax.set_ylabel("Frecuencia")
            ax.legend()     # Muestra la leyenda para identificar las líneas
            ax.grid(True, alpha=0.3)
            
            if USAR_LOG_SCALE:
                ax.set_yscale('log')

            # Guardar
            nombre_salida = f"{label}.png"
            ruta_salida = os.path.join(CARPETA_SALIDA, nombre_salida)
            
            plt.tight_layout()
            plt.savefig(ruta_salida, dpi=120)
            plt.close(fig)
            
            print(f"Guardado: {nombre_salida}")
            count += 1

        except Exception as e:
            print(f"Error en {nombre_archivo}: {e}")

    print(f"\n--- Listo. {count} gráficos guardados en '{CARPETA_SALIDA}' ---")

if __name__ == "__main__":
    procesar_grafico_unico()