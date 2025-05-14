import numpy as np
import matplotlib.pyplot as plt

# Helper function to read data from files
def read_columns(filename, usecols, delimiter=None, skiprows=0):
    data = np.loadtxt(filename, usecols=usecols, delimiter=delimiter, skiprows=skiprows)
    return data

# Read magnetizaciones.txt: columns: T, m_sup, m_inf, paso, m_avg
mag_data = read_columns('magnetizaciones.txt', usecols=(0,1,2,3,4))
T_mag = mag_data[:,0]
m_sup = mag_data[:,1]
m_inf = mag_data[:,2]
paso_mag = mag_data[:,3]
m_avg = mag_data[:,4]

# Read energias.txt: columns: T, E, paso
ener_data = read_columns('energias.txt', usecols=(0,1,2))
T_ener = ener_data[:,0]
E = ener_data[:,1]
paso_ener = ener_data[:,2]

# Lee histograma.txt (omite la primera línea)
hist_data = np.loadtxt('histograma.txt', skiprows=1)
spins_plus = hist_data[:,0]
cuentas = hist_data[:,1]

# Read chi.txt: columns: T, paso, chi
chi_data = read_columns('chi.txt', usecols=(0,1,2))
T_chi = chi_data[:,0]
paso_chi = chi_data[:,1]
chi = chi_data[:,2]

# Read cv.txt: columns: T, paso, cv
cv_data = read_columns('cv.txt', usecols=(0,1,2))
T_cv = cv_data[:,0]
paso_cv = cv_data[:,1]
cv = cv_data[:,2]

# Assume temperature is constant for all files (from first value)
main_T = T_mag[0]

fig, axs = plt.subplots(3, 2, figsize=(14, 12))
fig.suptitle(f'Ising Kawasaki Simulation Results at T = {main_T}', fontsize=16)

# Magnetization subplot
axs[0,0].plot(paso_mag, m_sup, label='Superior')
axs[0,0].plot(paso_mag, m_inf, label='Inferior')
axs[0,0].plot(paso_mag, m_avg, label='Media')
axs[0,0].set_title('Magnetización')
axs[0,0].set_xlabel('Paso')
axs[0,0].set_ylabel('Magnetización')
axs[0,0].legend()
axs[0,0].grid(True)

# Energy subplot
axs[0,1].plot(paso_ener, E, color='tab:red')
axs[0,1].set_title('Energía Total')
axs[0,1].set_xlabel('Paso')
axs[0,1].set_ylabel('Energía')
axs[0,1].grid(True)

# Histograma de barras de spins +1 por fila (red final)
axs[1,0].bar(spins_plus, cuentas, color='tab:green')
axs[1,0].set_title('Histograma de spins +1 por fila (red final)')
axs[1,0].set_xlabel('Número de spins +1 en una fila')
axs[1,0].set_ylabel('Número de filas')
axs[1,0].grid(True, axis='y')

# Susceptibility subplot
axs[1,1].plot(paso_chi, chi, color='tab:purple')
axs[1,1].set_title('Susceptibilidad Magnética')
axs[1,1].set_xlabel('Paso')
axs[1,1].set_ylabel('Susceptibilidad')
axs[1,1].grid(True)

# Specific heat subplot
axs[2,0].plot(paso_cv, cv, color='tab:orange')
axs[2,0].set_title('Calor Específico')
axs[2,0].set_xlabel('Paso')
axs[2,0].set_ylabel('Calor Específico')
axs[2,0].grid(True)

# Hide the empty subplot (bottom right)
axs[2,1].axis('off')

plt.tight_layout(rect=[0, 0.03, 1, 0.96])
plt.show()