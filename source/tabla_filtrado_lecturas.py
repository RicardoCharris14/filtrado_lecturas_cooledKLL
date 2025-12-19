import pandas as pd
import matplotlib.pyplot as plt
import os

# --- CONFIGURACIÓN ---
ARCHIVO_ENTRADA = 'data/filtrado/resultados_filtro_genomas.csv' 
# Nombre de la imagen de salida
ARCHIVO_SALIDA = 'data/filtrado/tabla_resultados_filtro.png'

# Mapeo para renombrar las columnas del CSV a algo más legible en la imagen
# Clave: Nombre en el CSV -> Valor: Nombre en la Tabla final
NOMBRES_COLUMNAS = {
    'k': 'K',
    'lower_quantile': 'Low Q',
    'upper_quantile': 'Up Q',
    'lower_bound': 'Low Bound',
    'upper_bound': 'Up Bound',
    'elements': 'Elements',
    'unique_elim_e': 'Unique Elim.',
    'elim_e': 'Total Elim.'
}

def generar_tabla_bounds():
    # 1. Crear datos de ejemplo si no existe el archivo (para que el script funcione al probarlo)
    if not os.path.exists(ARCHIVO_ENTRADA):
        print(f"{ARCHIVO_ENTRADA}, no existe.")
        exit(1)

    # 2. Leer el CSV
    try:
        df = pd.read_csv(ARCHIVO_ENTRADA)
    except Exception as e:
        print(f"Error al leer el CSV: {e}")
        return

    # 3. Formateo de datos
    # Formatear números grandes con separador de miles (ej: 1,000,000)
    cols_numericas = ['lower_bound', 'upper_bound', 'elements', 'unique_elim_e', 'elim_e']
    for col in cols_numericas:
        if col in df.columns:
            df[col] = df[col].apply(lambda x: f"{int(x):,}")

    # Renombrar columnas
    df_visual = df.rename(columns=NOMBRES_COLUMNAS)

    # 4. Generar la Tabla con Matplotlib
    # Ajustar tamaño de figura según filas y columnas
    fig_width = len(df_visual.columns) * 1.5
    fig_height = len(df_visual) * 0.5 + 1.5
    
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    
    # Ocultar ejes
    ax.axis('off')
    ax.axis('tight')

    # Crear tabla
    tabla = ax.table(
        cellText=df_visual.values,
        colLabels=df_visual.columns,
        cellLoc='center',
        loc='center'
    )

    # 5. Estilizar la tabla
    tabla.auto_set_font_size(False)
    tabla.set_fontsize(10)
    tabla.scale(1.0, 1.8) # Escalar celdas (Ancho, Alto)

    for (row, col), cell in tabla.get_celld().items():
        cell.set_edgecolor('white') # Bordes blancos para limpieza
        cell.set_linewidth(1)
        
        if row == 0:
            # Cabecera
            cell.set_text_props(weight='bold', color='white')
            cell.set_facecolor('#40466e') # Azul oscuro
            cell.set_height(0.1)
        else:
            # Filas de datos
            cell.set_height(0.08)
            if row % 2 == 0:
                cell.set_facecolor('#f2f2f2') # Gris claro alternado
            else:
                cell.set_facecolor('#ffffff') # Blanco
    
    # Título
    plt.title("Filtrado k-mers presentes en 50 genomas de Escherichia Coli", weight='bold', fontsize=14, pad=10)

    # 6. Guardar
    plt.savefig(ARCHIVO_SALIDA, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Tabla guardada exitosamente como: {ARCHIVO_SALIDA}")
    
    # Mostrar vista previa en texto
    print("\n--- Vista previa de los datos ---")
    print(df_visual.to_string(index=False))

if __name__ == "__main__":
    generar_tabla_bounds()