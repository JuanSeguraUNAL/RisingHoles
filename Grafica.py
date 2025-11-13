import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def cargar_trayectorias(nombre_archivo):
    trayectorias = []
    with open(nombre_archivo, 'r') as f:
        for linea in f:
            puntos = []
            for tripleta in linea.strip().replace(',', ' ').split():
                # Esperamos múltiplos de 3 valores
                continue
            # Alternativa robusta:
            coords = [float(x) for x in linea.replace(',', ' ').split()]
            puntos = np.array(coords).reshape(-1, 3)
            trayectorias.append(puntos)
    return trayectorias

def graficar_trayectorias(trayectorias, guardar=False, archivo_salida="trayectorias.png"):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    for puntos in trayectorias:
        xs, ys, zs = puntos[:,0], puntos[:,1], puntos[:,2]
        ax.plot(xs, ys, zs, lw=1.5, alpha=0.7)
        ax.scatter(xs[-1], ys[-1], zs[-1], color='red', s=10)  # punto final de cada burbuja

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Trayectorias de burbujas en la malla')
    ax.grid(True)
    plt.tight_layout()

    if guardar:
        plt.savefig(archivo_salida, dpi=300)
    plt.show()

if __name__ == "__main__":
    archivo = "trayectoria_final.dat"  # Ajusta el nombre si tu simulador guarda otro
    trayectorias = cargar_trayectorias(archivo)
    print(f"Se cargaron {len(trayectorias)} trayectorias.")
    graficar_trayectorias(trayectorias)
