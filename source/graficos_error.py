import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
# Importamos la librería para ajustar texto (asegúrate de instalarla: pip install adjustText)
try:
    from adjustText import adjust_text
    HAVE_ADJUST_TEXT = True
except ImportError:
    HAVE_ADJUST_TEXT = False

# --- CONFIGURACIÓN ---
SKETCH_SETTINGS = 'NB_100_BC_10_CS_100' 

# Usamos os.path.join para evitar problemas de / o \ en Windows
CARPETA_DATOS = os.path.join('data', 'frequency_distribution', SKETCH_SETTINGS)
CARPETA_SALIDA = os.path.join('data', 'images', SKETCH_SETTINGS)

print(f"--- Iniciando análisis ---")

# Crear carpeta de salida si no existe
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
        vector_mem = df_mem['vector_memory'].iloc[0]
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
            'vector_memory': vector_mem,
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

def plot_graph_dual(y_column_error, title, y_label_error, ax):
    """
    Grafica Error (Eje Izq) vs Memoria (Eje Der) con etiquetas inteligentes.
    """
    
    color_error = 'tab:blue'
    ax.set_xlabel('Sketch Memory (Bytes)')
    ax.set_ylabel(y_label_error, color=color_error, fontweight='bold')
    
    # Graficar la línea del Error
    line1 = ax.plot(df_results['sketch_memory'], df_results[y_column_error], 
                    marker='o', color=color_error, label='Sketch Error', linewidth=2)
    ax.tick_params(axis='y', labelcolor=color_error)
    ax.grid(True, alpha=0.3)

    # --- MAGIA PARA ETIQUETAS (NO SOLAPAMIENTO) ---
    texts = []
    for i in range(len(df_results)):
        row = df_results.iloc[i]
        # Usamos row['label'] si existe, o generamos una etiqueta genérica
        label_text = row['label'] if 'label' in row else f"{i}"
        
        t = ax.text(row['sketch_memory'], row[y_column_error], label_text, 
                    fontsize=9, color=color_error, fontweight='bold')
        texts.append(t)
    
    # adjust_text repele las etiquetas para que no se pisen
    if HAVE_ADJUST_TEXT:
        adjust_text(texts, ax=ax, 
                    arrowprops=dict(arrowstyle='-', color='gray', lw=0.5),
                    force_points=0.2, force_text=0.2, expand_points=(1.2, 1.2))

    # ---------------------------------------------------------
    # 2. Eje Derecho: MEMORIA (Vector y Comprimido)
    # ---------------------------------------------------------
    ax2 = ax.twinx()  
    color_vec = 'tab:red'
    color_comp = 'tab:green'
    ax2.set_ylabel('Memory Size (Bytes)', color='black', fontweight='bold')

    # Graficar Memoria Vector (x)
    line2 = ax2.plot(df_results['sketch_memory'], df_results['vector_memory'], 
                     marker='x', linestyle='--', color=color_vec, label='Vector Mem', alpha=0.7)
    
    # Graficar Memoria Comprimida (^)
    line3 = ax2.plot(df_results['sketch_memory'], df_results['compressed_memory'], 
                     marker='^', linestyle=':', color=color_comp, label='Compressed Mem', alpha=0.7)
    
    ax2.tick_params(axis='y', labelcolor='black')

    # ---------------------------------------------------------
    # 3. Leyenda Combinada
    # ---------------------------------------------------------
    lines = line1 + line2 + line3
    labels = [l.get_label() for l in lines]
    ax.legend(lines, labels, loc='best', fontsize='small', framealpha=0.9)

    ax.set_title(title)

# --- GENERACIÓN DE GRÁFICOS ---
fig, axs = plt.subplots(2, 2, figsize=(16, 12))
plt.subplots_adjust(hspace=0.3, wspace=0.35) 

# Gráficos
plot_graph_dual('mre_quantile', '1. Quantile: Relative Error vs Memory', 'MRE', axs[0, 0])
plot_graph_dual('mae_quantile', '2. Quantile: Absolute Error vs Memory', 'MAE', axs[0, 1])
plot_graph_dual('mre_rank', '3. Rank: Relative Error vs Memory', 'MRE', axs[1, 0])
plot_graph_dual('mae_rank', '4. Rank: Absolute Error vs Memory', 'MAE', axs[1, 1])

print("Generando gráficos actualizados...")

# Corrección segura de ruta para guardar
archivo_salida = os.path.join(CARPETA_SALIDA, 'memoryVSError.jpg')
plt.savefig(archivo_salida)
print(f"Gráfico guardado en: {archivo_salida}")

plt.show()