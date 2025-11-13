import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

datos = np.loadtxt("resultados_simulacion/estado_completo_paso_15.txt", skiprows=3)

print(f"Forma de los datos: {datos.shape}")

x, y, z, masa, temp, ocup = datos.T


fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
scatter = ax.scatter(x, y, z, c=masa, cmap='viridis', s=10, alpha=0.6)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Distribución de Masa 3D')
plt.colorbar(scatter, label='Masa')
plt.show()