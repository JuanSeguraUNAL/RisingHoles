import numpy as np
import plotly.graph_objects as go

def cargar_datos(archivo):
    # Carga del archivo (omitimos cabecera de 3 líneas)
    # altura radio angulo valor_influencia x y z
    datos = np.loadtxt(archivo, skiprows=4)
    h, r, theta, influencia, x, y, z, material = datos.T

    return x, y, z, influencia

def graficar_plotly(archivo):
    x, y, z, influencia = cargar_datos(archivo)

    fig = go.Figure(data=[
        go.Scatter3d(
            x=x, y=y, z=z,
            mode='markers',
            marker=dict(
                size=4,
                color=influencia,
                colorscale='Inferno',   # Escala de color
                colorbar=dict(title='Influencia'),
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
        title='Campo de influencia de burbujas en 3D',
        template='plotly_dark'
    )

    fig.show()

graficar_plotly("resultados_simulacion/influencia/influencia_paso_500.txt")