import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

def cargar_datos(archivo):
    datos = np.loadtxt(archivo, skiprows=4)
    h, r, theta, influencia, x, y, z, material = datos.T
    return r, x, y, z, influencia, material

def main():
    r, x, y, z, influencia, material = cargar_datos('resultados_simulacion/influencia/influencia_paso_500.txt')

    # --- Seleccionar capa superior ---
    z_max = np.max(z)

    mask = (z == z_max) 

    r_top = r[mask]
    x_top = x[mask] / np.max(r_top)
    y_top = y[mask] / np.max(r_top)
    influencia_top = influencia[mask]

    # --- Crear malla regular para imshow ---
    # Resolución de la malla (ajustable)
    grid_res = 200

    xi = np.linspace(x_top.min(), x_top.max(), grid_res)
    yi = np.linspace(y_top.min(), y_top.max(), grid_res)

    Xi, Yi = np.meshgrid(xi, yi)

    # Interpolar puntos a la malla
    Zi = griddata((x_top, y_top), influencia_top, (Xi, Yi), method='cubic')

    # --- Graficar con imshow ---
    plt.figure(figsize=(7, 6))

    plt.imshow(
        Zi,
        extent=(xi.min(), xi.max(), yi.min(), yi.max()),
        origin="lower",
        aspect="equal",
        cmap="plasma"
    )
    plt.colorbar(label="Influencia")

    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title(f"Scalar Bubble Field (z = {z_max})")
    plt.tight_layout()
    plt.savefig('Chimneys.png', dpi=300)
    plt.show()

    r_unique = np.unique(r_top)
    infl_mean = []

    for R in r_unique:
        infl_mean.append(np.mean(influencia_top[r_top == R]))

    fig, ax = plt.subplots(figsize=(10,6))
    ax.semilogy(r_unique / np.max(r_unique), infl_mean, 'k-', lw=2)
    ax.set_xlabel('Normalized radius', fontsize=14)
    ax.set_ylabel('Influence', fontsize=14)
    ax.set_xlim(0, 0.8)
    plt.show()

main()
