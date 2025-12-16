import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
import re

# --- CONFIGURACIÓN ---
CARPETA_EXPERIMENTOS = 'data/experiments'
CARPETA_SALIDA = 'data/images' 

# Suavizado (Media Móvil)
VENTANA_SUAVIZADO = 20

ESTRUCTURAS_INFO = {
    'plain':      {'label': 'Plain Vector',      'color': 'tab:blue',   'marker': 'o'},
    'compressed': {'label': 'Compressed Vector', 'color': 'tab:green',  'marker': 's'},
    'sketch':     {'label': 'Sketch',            'color': 'tab:red',    'marker': '^'}
}

def parsear_nombre_archivo(nombre_archivo):
    nombre_lower = nombre_archivo.lower()
    
    # 1. Identificar Estructura
    estructura_detectada = None
    for clave in ESTRUCTURAS_INFO.keys():
        if clave in nombre_lower:
            estructura_detectada = clave
            break
    
    if not estructura_detectada: return None, None

    # 2. Identificar K
    numeros = re.findall(r'\d+', nombre_archivo)
    if not numeros: return None, None
    
    try:
        k = int(numeros[-1])
    except:
        return None, None

    return estructura_detectada, k

def graficar_consultas_final():
    print(f"--- Escaneando carpeta: {CARPETA_EXPERIMENTOS} ---")
    print(f"--- Suavizado: {VENTANA_SUAVIZADO} ---")
    
    patron = os.path.join(CARPETA_EXPERIMENTOS, 'consultas*.csv')
    archivos = glob.glob(patron)
    
    if not archivos:
        print("[ERROR] No se encontraron archivos.")
        return

    datos_agrupados = {}
    
    for ruta in archivos:
        nombre_archivo = os.path.basename(ruta)
        struct, k = parsear_nombre_archivo(nombre_archivo)
        
        if struct and k is not None:
            if k not in datos_agrupados:
                datos_agrupados[k] = {}
            try:
                df = pd.read_csv(ruta, skipinitialspace=True)
                
                # --- SUAVIZADO ---
                if VENTANA_SUAVIZADO > 1:
                    cols = ['quantile_t_mean', 'quantile_t_stdev', 'rank_t_mean', 'rank_t_stdev']
                    cols_existentes = [c for c in cols if c in df.columns]
                    df[cols_existentes] = df[cols_existentes].rolling(
                        window=VENTANA_SUAVIZADO, min_periods=1, center=True
                    ).mean()

                datos_agrupados[k][struct] = df
            except Exception as e:
                print(f"  [ERROR] {nombre_archivo}: {e}")

    lista_k = sorted(datos_agrupados.keys())
    if not lista_k: return

    carpeta_salida = os.path.join(CARPETA_SALIDA, 'graficos_consultas')
    if not os.path.exists(carpeta_salida):
        os.makedirs(carpeta_salida)

    # Función auxiliar de graficado (Sombra)
    def plot_con_sombra(ax, x, y, err, estilo, linestyle='-'):
        ax.plot(
            x, y, 
            label=estilo['label'], 
            color=estilo['color'], 
            linewidth=1.5,
            linestyle=linestyle
        )
        ax.fill_between(
            x, y - err, y + err, 
            color=estilo['color'], alpha=0.2, edgecolor=None
        )

    # Función para dibujar en un par de ejes (Quantile y Rank)
    def dibujar_fila(ax_q, ax_r, grupo_datos, excluir_plain=False):
        # --- QUANTILE ---
        hay_datos_q = False
        for struct, df in grupo_datos.items():
            if excluir_plain and struct == 'plain': continue
            if 'quantile_t_mean' not in df.columns: continue
            
            estilo = ESTRUCTURAS_INFO[struct]
            t_mean = df['quantile_t_mean']
            t_std = df['quantile_t_stdev']
            
            plot_con_sombra(ax_q, df['quantile'], t_mean, t_std, estilo)
            hay_datos_q = True

        ax_q.set_xlabel("Quantile")
        ax_q.set_ylabel("Tiempo (ns)")
        if hay_datos_q:
            ax_q.legend()
            ax_q.grid(True, alpha=0.3)

        # --- RANK ---
        hay_datos_r = False
        for struct, df in grupo_datos.items():
            if excluir_plain and struct == 'plain': continue
            if 'rank_t_mean' not in df.columns: continue
            
            estilo = ESTRUCTURAS_INFO[struct]
            t_mean = df['rank_t_mean']
            t_std = df['rank_t_stdev']
            
            plot_con_sombra(ax_r, df['quantile'], t_mean, t_std, estilo, linestyle='--')
            hay_datos_r = True

        ax_r.set_xlabel("Quantile")
        ax_r.set_ylabel("Tiempo (ns)")
        if hay_datos_r:
            ax_r.legend()
            ax_r.grid(True, alpha=0.3)

    # --- BUCLE PRINCIPAL POR K ---
    for k in lista_k:
        grupo = datos_agrupados[k]
        print(f"Generando figura combinada K={k}...")
        
        # Creamos una figura de 2 filas x 2 columnas
        fig, axs = plt.subplots(2, 2, figsize=(16, 12))
        
        # Título General
        suffix_title = f" - Suavizado ({VENTANA_SUAVIZADO})" if VENTANA_SUAVIZADO > 1 else ""
        fig.suptitle(f'Comparativa Completa Consultas (K = {k}){suffix_title}', fontsize=16, fontweight='bold')

        # === FILA 1: TODOS (Con Plain) ===
        axs[0, 0].set_title("Quantile vs Tiempo", fontweight='bold')
        axs[0, 1].set_title("Rank vs Tiempo", fontweight='bold')
        dibujar_fila(axs[0, 0], axs[0, 1], grupo, excluir_plain=False)

        # === FILA 2: ZOOM (Sin Plain) ===
        axs[1, 0].set_title("Quantile vs Tiempo", fontweight='bold')
        axs[1, 1].set_title("Rank vs Tiempo", fontweight='bold')
        dibujar_fila(axs[1, 0], axs[1, 1], grupo, excluir_plain=True)

        plt.tight_layout()
        # Ajustamos el espaciado superior para que no tape el título principal
        plt.subplots_adjust(top=0.92)
        
        nombre_salida = f"Consultas_K{k}.png"
        ruta_final = os.path.join(carpeta_salida, nombre_salida)
        
        plt.savefig(ruta_final, dpi=150)
        plt.close(fig)
        print(f"  -> Guardado: {ruta_final}")

    print(f"\n--- Listo. Gráficos en: {carpeta_salida} ---")

if __name__ == "__main__":
    graficar_consultas_final()