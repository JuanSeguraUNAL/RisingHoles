import os
import numpy as np
import seaborn as sns

# Diseño de los gráficos
sns.set()
sns.set_context("paper")
sns.set_palette("colorblind")

def cargar_datos(archivo):
    # Carga del archivo (omitimos cabecera de 3 líneas)
    # altura radio angulo valor_temperatura x y z
    datos = np.loadtxt(archivo, skiprows=4)
    h, r, theta, temperatura, x, y, z = datos.T

    return h, r, theta, x, y, z, temperatura

def 

