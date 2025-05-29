import numpy as np
import matplotlib.pyplot as plt

# Lista de tamaños de red disponibles
Ns = [32, 64, 128]

# Helper para leer dos columnas
def read_two_columns(filename):
    return np.loadtxt(filename, unpack=True)

fig, axs = plt.subplots(3, 3, figsize=(18, 12))
fig.suptitle('Ising Kawasaki Simulation Results (comparativa tamaños de red)', fontsize=16)

# --- Calor específico (Cv) ---
for N in Ns:
    try:
        T_cv, cv = read_two_columns(f'cv{N}.txt')
        axs[2,0].plot(T_cv, cv, label=f'N={N}')
    except Exception as e:
        print(f'No se pudo leer cv{N}.txt: {e}')
axs[2,0].set_title('Calor Específico')
axs[2,0].set_xlabel('Temperatura')
axs[2,0].set_ylabel('Cv')
axs[2,0].legend()
axs[2,0].grid(True)

# --- Susceptibilidad (chi) ---
for N in Ns:
    try:
        T_chi, chi = read_two_columns(f'chi{N}.txt')
        axs[1,2].plot(T_chi, chi, label=f'N={N}')
    except Exception as e:
        print(f'No se pudo leer chi{N}.txt: {e}')
axs[1,2].set_title('Susceptibilidad Magnética')
axs[1,2].set_xlabel('Temperatura')
axs[1,2].set_ylabel('Chi')
axs[1,2].legend()
axs[1,2].grid(True)

# --- Densidad ---
for N in Ns:
    try:
        T_dens, dens = read_two_columns(f'densidades{N}.txt')
        axs[1,0].plot(T_dens, dens, label=f'N={N}')
    except Exception as e:
        print(f'No se pudo leer densidades{N}.txt: {e}')
axs[1,0].set_title('Densidad promedio')
axs[1,0].set_xlabel('Temperatura')
axs[1,0].set_ylabel('Densidad')
axs[1,0].legend()
axs[1,0].grid(True)

# --- Energía ---
for N in Ns:
    try:
        T_ener, ener = read_two_columns(f'energias{N}.txt')
        axs[2,1].plot(T_ener, ener, label=f'N={N}')
    except Exception as e:
        print(f'No se pudo leer energias{N}.txt: {e}')
axs[2,1].set_title('Energía promedio')
axs[2,1].set_xlabel('Temperatura')
axs[2,1].set_ylabel('Energía')
axs[2,1].legend()
axs[2,1].grid(True)

# --- Magnetización: una subgráfica para cada N ---
for idx, N in enumerate(Ns):
    try:
        T_mag, m_sup, m_inf, m_med = np.loadtxt(f'magnetizaciones{N}.txt', unpack=True)
        axs[0,idx].plot(T_mag, m_sup, label='Sup')
        axs[0,idx].plot(T_mag, m_inf, label='Inf')
        axs[0,idx].plot(T_mag, m_med, label='Med')
        axs[0,idx].set_title(f'Magnetización N={N}')
        axs[0,idx].set_xlabel('Temperatura')
        axs[0,idx].set_ylabel('Magnetización')
        axs[0,idx].legend()
        axs[0,idx].grid(True)
    except Exception as e:
        print(f'No se pudo leer magnetizaciones{N}.txt: {e}')
        axs[0,idx].set_visible(False)

# --- Subgráficas vacías ---
axs[1,1].axis('off')
axs[2,2].axis('off')

plt.tight_layout(rect=[0, 0.03, 1, 0.96])
plt.show()