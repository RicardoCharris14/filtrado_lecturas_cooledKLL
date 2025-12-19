import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import numpy as np

# --- CONFIGURACIÓN ---
SKETCH_SETTINGS = 'NB_100_BC_10_CS_100' # !!!Modificar para graficar los datos requeridos!!!

CARPETA_DATOS = 'data/frequency_distribution/'+SKETCH_SETTINGS        # Dónde están los CSV
CARPETA_SALIDA = 'data/images/'+SKETCH_SETTINGS # Nombre de la carpeta donde se guardarán las imágenes

# CONFIGURACIÓN DE MEDIA MÓVIL
USAR_MEDIA_MOVIL = True
TAMAÑO_VENTANA = 5

print(f"--- Iniciando generación de gráficos ---\nLeeyendo datos de: {CARPETA_DATOS}")

# Crear carpeta de salida si no existe
if not os.path.exists(CARPETA_SALIDA):
    os.makedirs(CARPETA_SALIDA)
    print(f"Carpeta creada: {CARPETA_SALIDA}")

patron = os.path.join(CARPETA_DATOS, '*mers_distribution.csv')
dist_files = glob.glob(patron)
dist_files.sort()

if not dist_files:
    print("[ERROR] No se encontraron archivos CSV.")
    exit()

# Función de graficado reutilizable
def plot_dual_axis(ax, x, y1, y2, title, y1_label, y2_label, label1, label2):
    # Eje Izquierdo
    color1 = 'tab:blue'
    ax.set_xlabel('Quantile (0.0 - 1.0)')
    ax.set_ylabel(y1_label, color=color1, fontweight='bold')
    line1 = ax.plot(x, y1, color=color1, label=label1, linewidth=1.5)
    ax.tick_params(axis='y', labelcolor=color1)
    ax.grid(True, alpha=0.3)

    # Eje Derecho
    ax_right = ax.twinx()
    color2 = 'tab:orange'
    ax_right.set_ylabel(y2_label, color=color2, fontweight='bold')
    line2 = ax_right.plot(x, y2, color=color2, label=label2, linewidth=1.5, linestyle='--')
    ax_right.tick_params(axis='y', labelcolor=color2)

    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax.legend(lines, labels, loc='upper center', fontsize='small', framealpha=0.9)
    ax.set_title(title)

for dist_file in dist_files:
    try:
        nombre_archivo = os.path.basename(dist_file)
        label = nombre_archivo.replace('_distribution.csv', '')
        
        df = pd.read_csv(dist_file)
        if df.empty: continue

        # --- NUEVA LÓGICA: LEER MEMORY CSV ---
        # Construimos la ruta del archivo de memoria reemplazando el sufijo
        memory_file = dist_file.replace('distribution.csv', 'memory.csv')
        info_unique = ""

        if os.path.exists(memory_file):
            try:
                df_mem = pd.read_csv(memory_file)
                # Verificamos que exista la columna y el archivo no esté vacío
                if 'unique_elements' in df_mem.columns and not df_mem.empty:
                    val_unique = df_mem.iloc[0]['unique_elements']
                    val_total = df_mem.iloc[0]['elements']
                    info_unique = f" - Elementos total/unicos: {val_total}/{val_unique}"
                else:
                    print(f"[AVISO] Columna 'unique_elements' no encontrada en {os.path.basename(memory_file)}")
            except Exception as e:
                print(f"[ERROR] Al leer memoria {os.path.basename(memory_file)}: {e}")
        else:
            print(f"[AVISO] No se encontró archivo de memoria para {label}")
        # -------------------------------------

        # --- CÁLCULOS ---
        df['abs_err_quantile'] = (df['real_quantile'] - df['estimated_quantile']).abs()
        df['abs_err_rank'] = (df['real_rank'] - df['estimated_rank']).abs()

        df['rel_err_quantile'] = np.where(df['real_quantile'] != 0, 
                                          df['abs_err_quantile'] / df['real_quantile'], 0)
        df['rel_err_rank'] = np.where(df['real_rank'] != 0, 
                                      df['abs_err_rank'] / df['real_rank'], 0)

        # --- SUAVIZADO ---
        if USAR_MEDIA_MOVIL:
            cols_to_smooth = ['abs_err_quantile', 'abs_err_rank', 'rel_err_quantile', 'rel_err_rank']
            for col in cols_to_smooth:
                df[col] = df[col].rolling(window=TAMAÑO_VENTANA, min_periods=1, center=True).mean()
            titulo_extra = f"(Suavizado: {TAMAÑO_VENTANA})"
        else:
            titulo_extra = ""

        # --- GRAFICACIÓN ---
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        
        # Título actualizado con la información de unique_elements
        fig.suptitle(f'Análisis: {label} {titulo_extra}{info_unique}', fontsize=16, fontweight='bold')

        # Gráfico 1: Relativo
        plot_dual_axis(ax1, df['quantile'], df['rel_err_quantile'], df['rel_err_rank'],
                       "Error Relativo vs Quantile", "Rel. Error: Quantile", "Rel. Error: Rank",
                       "Quantile Err (Rel)", "Rank Err (Rel)")

        # Gráfico 2: Absoluto
        plot_dual_axis(ax2, df['quantile'], df['abs_err_quantile'], df['abs_err_rank'],
                       "Error Absoluto vs Quantile", "Abs. Error: Quantile", "Abs. Error: Rank",
                       "Quantile Err (Abs)", "Rank Err (Abs)")

        plt.tight_layout()
        
        # --- GUARDAR ---
        nombre_imagen = f"{label}_analisis.png"
        ruta_completa = os.path.join(CARPETA_SALIDA, nombre_imagen)
        
        plt.savefig(ruta_completa, dpi=150) 
        print(f"Guardado: {ruta_completa}")
        
        plt.close(fig)

    except Exception as e:
        print(f"Error procesando {label}: {e}")

print(f"\n--- Proceso completado. Revisa la carpeta '{CARPETA_SALIDA}' ---")