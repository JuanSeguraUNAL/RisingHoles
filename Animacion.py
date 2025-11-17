import os
import numpy as np
import plotly.graph_objects as go

# === CONFIGURACIÓN ===
carpeta = "resultados_simulacion/influencia"
archivos = sorted([
    os.path.join(carpeta, f)
    for f in os.listdir(carpeta)
    if f.startswith("influencia_paso_") and f.endswith(".txt")
])

def cargar_datos_completos(archivo):
    """Carga el archivo COMPLETO y retorna coordenadas + valor de influencia."""
    try:
        print(f"  Cargando {os.path.basename(archivo)}...")
        
        # Leer TODAS las líneas (sin límite)
        with open(archivo, 'r') as f:
            lines = []
            for i, line in enumerate(f):
                if i < 4:  # Saltar encabezados
                    continue
                if line.strip() and not line.startswith('#'):
                    lines.append(line)
        
        print(f"    {len(lines)} líneas de datos encontradas")
        
        # Procesar datos
        datos = []
        for line in lines:
            parts = line.strip().split()
            if len(parts) >= 7:
                datos.append([float(p) for p in parts[:7]])
        
        datos = np.array(datos)
        if len(datos) == 0:
            print(f"    ¡ADVERTENCIA: No hay datos en {archivo}!")
            return np.array([]), np.array([]), np.array([]), np.array([])
            
        h, r, theta, influencia, x, y, z = datos.T
        
        print(f"    {len(x)} puntos cargados")
        print(f"    Rango Z: {np.min(z):.1f} a {np.max(z):.1f}")
        print(f"    Rango influencia: {np.min(influencia):.3f} a {np.max(influencia):.3f}")
        
        return x, y, z, influencia
        
    except Exception as e:
        print(f"ERROR cargando {archivo}: {e}")
        return np.array([]), np.array([]), np.array([]), np.array([])

# === Calcular rango global de influencia ===
print("Calculando rango global de influencia (cargando TODOS los datos)...")
inf_min_global = float('inf')
inf_max_global = float('-inf')
z_min_global = float('inf')
z_max_global = float('-inf')
x_range_global = 0
y_range_global = 0

# Usar solo algunos archivos para el cálculo del rango (para ahorrar tiempo)
archivos_calculo_rango = archivos[:5]  # Primeros 5 archivos para calcular rangos
print(f"Usando {len(archivos_calculo_rango)} archivos para cálculo de rangos")

for i, archivo in enumerate(archivos_calculo_rango):
    print(f"Procesando archivo {i+1}/{len(archivos_calculo_rango)} para rangos...")
    x, y, z, inf = cargar_datos_completos(archivo)
    if len(inf) > 0:
        inf_min_global = min(inf_min_global, np.min(inf))
        inf_max_global = max(inf_max_global, np.max(inf))
        z_min_global = min(z_min_global, np.min(z))
        z_max_global = max(z_max_global, np.max(z))
        x_range_global = max(x_range_global, np.max(np.abs(x)))
        y_range_global = max(y_range_global, np.max(np.abs(y)))

# Si no se encontraron datos, usar valores por defecto
if inf_min_global == float('inf'):
    inf_min_global = 0
    inf_max_global = 1
    z_min_global = 0
    z_max_global = 25
    x_range_global = 15
    y_range_global = 15

print(f"\nRango global de influencia: {inf_min_global:.3f} a {inf_max_global:.3f}")
print(f"Rango global Z: {z_min_global:.1f} a {z_max_global:.1f}")

# Asegurar que el eje Z sea más prominente
z_range = z_max_global - z_min_global
xy_range = max(x_range_global, y_range_global) * 2

# Ajustar el ratio para hacer el eje Z más grande
z_scale_factor = 2.0  # Puedes ajustar este valor
aspect_ratio = dict(x=1, y=1, z=z_scale_factor * (z_range / xy_range) if xy_range > 0 else 1)

print(f"Rango XY: ±{xy_range/2:.1f}")
print(f"Ratio de aspecto: {aspect_ratio}")

# === Cargar el primer frame COMPLETO ===
print(f"\nCargando primer frame COMPLETO: {os.path.basename(archivos[0])}")
x0, y0, z0, inf0 = cargar_datos_completos(archivos[0])

if len(x0) == 0:
    print("ERROR: No se pudieron cargar datos del primer archivo.")
    exit()

print(f"✓ Primer frame cargado: {len(x0)} puntos")

# Normalizar el tamaño de los marcadores (usando rango global)
tamano_min, tamano_max = 2, 15
if inf_max_global > inf_min_global:
    tamanos = tamano_min + (tamano_max - tamano_min) * (inf0 - inf_min_global) / (inf_max_global - inf_min_global)
else:
    tamanos = np.ones_like(inf0) * tamano_min

# === Crear figura base ===
fig = go.Figure()

# Añadir puntos con color y tamaño según influencia (usando escala fija)
scatter = go.Scatter3d(
    x=x0, y=y0, z=z0,
    mode='markers',
    marker=dict(
        size=tamanos,
        color=inf0,
        colorscale='Inferno',
        cmin=inf_min_global,  # ESCALA FIJA
        cmax=inf_max_global,  # ESCALA FIJA
        colorbar=dict(title="Influencia"),
        opacity=0.7,  # Un poco menos de opacidad para ver mejor la densidad
        showscale=True
    ),
    name="Trayectorias",
    hovertemplate="<b>Influencia:</b> %{marker.color:.2f}<br>" +
                  "X: %{x:.2f}<br>Y: %{y:.2f}<br>Z: %{z:.2f}<extra></extra>"
)

fig.add_trace(scatter)

# Configurar layout con eje Z más grande
fig.update_layout(
    title=f"Evolución del Campo de Influencia - {os.path.basename(archivos[0])}",
    scene=dict(
        xaxis_title="X (cm)",
        yaxis_title="Y (cm)", 
        zaxis_title="Z (cm)",
        aspectratio=aspect_ratio,  # Hace el eje Z más grande
        camera=dict(
            eye=dict(x=1.5, y=1.5, z=1.5 * z_scale_factor)  # Ajustar cámara para el Z más grande
        ),
        # Establecer rangos fijos para mejor visualización
        xaxis=dict(range=[-xy_range/2, xy_range/2]),
        yaxis=dict(range=[-xy_range/2, xy_range/2]),
        zaxis=dict(range=[z_min_global, z_max_global])
    ),
    margin=dict(l=0, r=0, b=0, t=40),
    updatemenus=[dict(
        type="buttons",
        showactive=False,
        buttons=[
            dict(label="▶️ Play", method="animate",
                 args=[None, {"frame": {"duration": 300, "redraw": True},  # Más tiempo por frame
                              "fromcurrent": True, "mode": "immediate"}]),
            dict(label="⏸️ Pause", method="animate",
                 args=[[None], {"frame": {"duration": 0, "redraw": False},
                                "mode": "immediate"}])
        ],
        x=0.1, y=0, xanchor="right", yanchor="bottom"
    )]
)

# === Crear los frames COMPLETOS ===
print(f"\nCreando {len(archivos)-1} frames de animación (esto puede tomar tiempo)...")
frames = []

for i, archivo in enumerate(archivos[1:]):  # Empezar desde el segundo archivo
    print(f"Creando frame {i+1}/{len(archivos)-1}: {os.path.basename(archivo)}")
    
    x, y, z, inf = cargar_datos_completos(archivo)
    
    if len(x) == 0:
        print(f"  ⚠️ Saltando frame {i+1} - sin datos")
        continue
        
    # Normalizar tamaños para este frame (usando rango global)
    if inf_max_global > inf_min_global:
        tamanos_frame = tamano_min + (tamano_max - tamano_min) * (inf - inf_min_global) / (inf_max_global - inf_min_global)
    else:
        tamanos_frame = np.ones_like(inf) * tamano_min
    
    frame = go.Frame(
        data=[go.Scatter3d(
            x=x, y=y, z=z,
            mode='markers',
            marker=dict(
                size=tamanos_frame,
                color=inf,
                colorscale='Inferno',
                cmin=inf_min_global,  # ESCALA FIJA
                cmax=inf_max_global,  # ESCALA FIJA
                opacity=0.7
            ),
            hovertemplate="<b>Influencia:</b> %{marker.color:.2f}<br>" +
                         "X: %{x:.2f}<br>Y: %{y:.2f}<br>Z: %{z:.2f}<extra></extra>"
        )],
        name=os.path.basename(archivo)
    )
    frames.append(frame)
    
    print(f"  ✓ Frame {i+1} creado con {len(x)} puntos")

fig.frames = frames

# Añadir slider
steps = []
for i, archivo in enumerate(archivos):
    step = dict(
        method="animate",
        args=[[f"frame_{i}"], {"mode": "immediate", "frame": {"duration": 300, "redraw": True}}],
        label=os.path.basename(archivo).replace('influencia_paso_', '').replace('.txt', '')
    )
    steps.append(step)

sliders = [dict(
    steps=steps,
    active=0,
    currentvalue={"prefix": "Paso: "},
    x=0.1, y=0,
    xanchor="left",
    len=0.9
)]

fig.update_layout(sliders=sliders)

print(f"\n🎉 Animación completada!")
print(f"📊 Estadísticas:")
print(f"   - {len(archivos)} archivos procesados")
print(f"   - {len(frames)} frames creados")
print(f"   - Escala de color fija: {inf_min_global:.2f} a {inf_max_global:.2f}")
print(f"   - Dimensión Z ampliada: factor {z_scale_factor}")
print(f"   - Puntos en primer frame: {len(x0)}")
print(f"\n📈 Mostrando animación...")

fig.show()