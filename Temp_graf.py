import numpy as np
import plotly.graph_objects as go

def cargar_datos(archivo):
    # Carga del archivo (omitimos cabecera de 3 líneas)
    # altura radio angulo valor_influencia x y z
    datos = np.loadtxt(archivo, skiprows=4)
    h, r, theta, temperatura, alpha, x, y, z, material = datos.T

    return x, y, z, temperatura

def graficar_plotly(archivo):
    x, y, z, temperatura= cargar_datos(archivo)

    fig = go.Figure(data=[
        go.Scatter3d(
            x=x, y=y, z=z,
            mode='markers',
            marker=dict(
                size=4,
                color=temperatura,
                colorscale='Inferno',   # Escala de color
                colorbar=dict(title='Temperature [°C]'),
                opacity=0.8
            )
        )
    ])

    fig.update_layout(
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z',
            aspectmode='cube'
        ),
        title='Temperature',
        template='plotly_dark'
    )

    fig.show()

graficar_plotly("resultados_simulacion/temperatura/temperatura_paso_500.txt")