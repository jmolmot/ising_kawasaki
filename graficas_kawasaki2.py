import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# Busca todos los N disponibles para un tipo de archivo
def buscar_N(archivo_base):
    archivos = glob.glob(f"{archivo_base}*.txt")
    Ns = []
    for f in archivos:
        nombre = os.path.splitext(os.path.basename(f))[0]
        n_str = nombre.replace(archivo_base, "")
        try:
            Ns.append(int(n_str))
        except ValueError:
            pass
    return sorted(Ns)

# Diccionario de magnitudes y archivos base
archivos = {
    "Mag. Superior": "promedio_magnetizacionsuperior",
    "Mag. Inferior": "promedio_magnetizacioninferior",
    "Densidad +": "promedio_densidadpositivo",
    "Densidad -": "promedio_densidadnegativo",
    "Calor específico": "filecv",
}
archivo_chi = "susceptibilidad"
archivo_energia = "energia_pmontecarlo"

# --- Subgráficas para magnitudes principales ---
fig, axs = plt.subplots(len(archivos), 1, figsize=(8, 12), sharex=True)
for idx, (titulo, base) in enumerate(archivos.items()):
    Ns = buscar_N(base)
    for N in Ns:
        try:
            y, x = np.loadtxt(f"{base}{N}.txt", unpack=True)
            axs[idx].plot(x, y, marker='o', label=f"N={N}")
        except Exception as e:
            print(f"No se pudo leer {base}{N}.txt: {e}")
    axs[idx].set_ylabel(titulo)
    axs[idx].legend()
    axs[idx].grid(True)
axs[-1].set_xlabel("Temperatura")
for ax in axs:
    ax.grid(True)
plt.tight_layout()
plt.suptitle("Magnitudes vs Temperatura para distintos N", y=1.02)
plt.show()

# --- Susceptibilidad magnética ---
plt.figure(figsize=(8,5))
Ns = buscar_N(archivo_chi)
for N in Ns:
    try:
        y, x = np.loadtxt(f"{archivo_chi}{N}.txt", unpack=True)
        plt.plot(x, y, marker='o', label=f"N={N}")
    except Exception as e:
        print(f"No se pudo leer {archivo_chi}{N}.txt: {e}")
plt.ylabel("Susceptibilidad magnética")
plt.xlabel("Temperatura")
plt.legend()
plt.title("Susceptibilidad magnética vs Temperatura")
plt.grid(True)
plt.tight_layout()
plt.show()