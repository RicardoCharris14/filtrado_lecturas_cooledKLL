import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

# Importamos la librería para ajustar texto
try:
    from adjustText import adjust_text
    HAVE_ADJUST_TEXT = True
except ImportError:
    HAVE_ADJUST_TEXT = False

# --- CONFIGURACIÓN ---
SKETCH_SETTINGS = 'NB_100_BC_10_CS_100' 

# Construcción de rutas segura para el sistema operativo
CARPETA_DATOS = os.path.join('data', 'frequency_distribution', SKETCH_SETTINGS)
CARPETA_SALIDA = os.path.join('data', 'images', SKETCH_SETTINGS)

print(f"--- Iniciando análisis (Solo Memoria Comprimida) ---")

# Crear carpeta de salida si no existe (exist_ok=True evita errores si ya existe)
if not os.path.exists(CARPETA_SALIDA):
    os.makedirs(CARPETA_SALIDA, exist_ok=True)
    print(f"Carpeta creada: {CARPETA_SALIDA}")

patron = os.path.join(CARPETA_DATOS, '*mers_distribution.csv')
dist_files = glob.glob(patron)

if not dist_files:
    print(f"[ERROR] No se encontraron archivos en: {CARPETA_DATOS}")
    exit()

data_points = []

for dist_file in dist_files:
    try:
        # Preparar nombres de archivo
        memory_file = dist_file.replace('distribution.csv', 'memory.csv')
        nombre_archivo = os.path.basename(dist_file)
        label = nombre_archivo.replace('_distribution.csv', '')

        if not os.path.exists(memory_file):
            continue

        # Leer CSVs
        df_dist = pd.read_csv(dist_file)
        df_mem = pd.read_csv(memory_file)

        if df_dist.empty or df_mem.empty:
            continue

        # Datos de memoria
        sketch_mem = df_mem['sketch_memory'].iloc[0]
        # vector_mem = df_mem['vector_memory'].iloc[0] # No se usa en este script
        compressed_mem = df_mem['compressed_vector_memory'].iloc[0]

        # --- Cálculos de Error ---
        # Quantile
        mae_quantile = (df_dist['real_quantile'] - df_dist['estimated_quantile']).abs().mean()
        mask_q = df_dist['real_quantile'] != 0
        mre_quantile = ((df_dist.loc[mask_q, 'real_quantile'] - df_dist.loc[mask_q, 'estimated_quantile']).abs() / df_dist.loc[mask_q, 'real_quantile']).mean() if mask_q.any() else 0

        # Rank
        mae_rank = (df_dist['real_rank'] - df_dist['estimated_rank']).abs().mean()
        mask_r = df_dist['real_rank'] != 0
        mre_rank = ((df_dist.loc[mask_r, 'real_rank'] - df_dist.loc[mask_r, 'estimated_rank']).abs() / df_dist.loc[mask_r, 'real_rank']).mean() if mask_r.any() else 0

        data_points.append({
            'label': label,
            'sketch_memory': sketch_mem,
            'compressed_memory': compressed_mem,
            'mae_quantile': mae_quantile,
            'mre_quantile': mre_quantile,
            'mae_rank': mae_rank,
            'mre_rank': mre_rank
        })

    except Exception as e:
        print(f"Error en {label}: {e}")

df_results = pd.DataFrame(data_points)
if df_results.empty:
    print("No hay datos para graficar.")
    exit()

# Ordenar por Sketch Memory
df_results = df_results.sort_values(by='sketch_memory')

# --- FUNCIÓN DE GRAFICADO ACTUALIZADA (Evita líneas) ---
def plot_graph_dual(y_column_error, title, y_label_error, ax):
    
    # 1. Eje Izquierdo: ERROR (Sketch)
    color_error = 'tab:blue'
    ax.set_xlabel('Sketch Memory (Bytes)')
    ax.set_ylabel(y_label_error, color=color_error, fontweight='bold')
    
    # Guardamos la línea en una variable (line1 es una lista de líneas)
    line1 = ax.plot(df_results['sketch_memory'], df_results[y_column_error], 
                    marker='o', color=color_error, label='Sketch Error', linewidth=2)
    ax.tick_params(axis='y', labelcolor=color_error)
    ax.grid(True, alpha=0.3)

    # 2. Eje Derecho: MEMORIA (Solo Compressed)
    ax2 = ax.twinx()
    color_comp = 'tab:green'
    ax2.set_ylabel('Compressed Memory (Bytes)', color=color_comp, fontweight='bold')

    # Guardamos la segunda línea
    line2 = ax2.plot(df_results['sketch_memory'], df_results['compressed_memory'], 
                     marker='^', linestyle=':', color=color_comp, label='Compressed Mem', alpha=0.7)
    
    ax2.tick_params(axis='y', labelcolor=color_comp)

    # Combinar leyendas
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax.legend(lines, labels, loc='upper left', fontsize='small')
    ax.set_title(title)

    # --- CORRECCIÓN DE ETIQUETAS ---
    texts = []
    for i in range(len(df_results)):
        row = df_results.iloc[i]
        # Creamos el texto
        t = ax.text(row['sketch_memory'], row[y_column_error], row['label'], 
                    fontsize=8, color=color_error, fontweight='bold')
        texts.append(t)
    
    if HAVE_ADJUST_TEXT:
        # AQUI ESTA EL CAMBIO CLAVE:
        # 1. add_objects: Le decimos que evite pisar las líneas (line1 + line2)
        # 2. expand_points: Aumentamos el "área" de los puntos para empujar el texto más lejos (factor 1.5 o 2)
        adjust_text(texts, 
                    ax=ax, 
                    add_objects=line1 + line2,  # ¡Evita las líneas!
                    expand_points=(1.5, 1.5),   # Empuja el texto más lejos del punto
                    arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

# --- GENERACIÓN DE GRÁFICOS ---
fig, axs = plt.subplots(2, 2, figsize=(16, 12))
plt.subplots_adjust(hspace=0.3, wspace=0.35) 

# Gráficos
plot_graph_dual('mre_quantile', '1. Quantile: Relative Error vs Memory', 'MRE', axs[0, 0])
plot_graph_dual('mae_quantile', '2. Quantile: Absolute Error vs Memory', 'MAE', axs[0, 1])
plot_graph_dual('mre_rank', '3. Rank: Relative Error vs Memory', 'MRE', axs[1, 0])
plot_graph_dual('mae_rank', '4. Rank: Absolute Error vs Memory', 'MAE', axs[1, 1])

print("Generando gráficos limpios (sin vector original)...")

# Guardado seguro
archivo_salida = os.path.join(CARPETA_SALIDA, 'memoryVSErrorZoomed.jpg')
plt.savefig(archivo_salida)
print(f"Gráfico guardado en: {archivo_salida}")

plt.show()