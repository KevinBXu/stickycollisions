import cupy as np
import math
import matplotlib.pyplot as plt
import scipy as sc
import importlib
import sys
sys.path.append('..')
import collision_cp
importlib.reload(collision_cp)
from collision_cp import generate_spectrum, wigner_smith_matrix, wigner_smith_matrix_single, thermal_time_delay, relative_MB_weights, S_matrix

# parameters
e_mass = 5.48579909e-4    # electron mass in amu
t0 = 315775               # a.u. to K
tau0 = 2.418884326509e-17 # a.u. to sec

# define properties of molecules
reduced_masses = {"NaKNaK" : 62 / 2, "RbKRb" : 87 * (87 + 40) / (87 + 87 + 40)}     # amu
van_der_wall_coeffs = {"NaKNaK" : 561070, "RbKRb" : 8000}      # a.u.

molecules = "RbCsRbCs"
mu = 110 / e_mass      # now in a.u.
c6 = 1.9e5         # a.u.

# define natural units
beta = (2 * mu * c6) ** (1/4)    # length 
E_beta = 1 / (2 * mu * beta ** 2)   # energy
tau_beta = 1 / E_beta    # time

num_samples = 5000

scale = 1

# properties of resonant spectrum
rrkm = 0.25e-3     # seconds
dos = rrkm / tau0 / 2 / np.pi * E_beta
mean_spacing = 1 / dos / scale

# thermal properties
temp = 1 / scale       # in terms of E6
num_resonances = 400        # should be odd

rrkm_c6 = 2 * np.pi / mean_spacing 

mean_x = 1
mean_coupling_strength = np.sqrt(mean_x * mean_spacing / np.pi ** 2)

energy_GOE_mat, W_mat, x_actual = generate_spectrum(num_resonances, mean_spacing, mean_coupling_strength)

from scipy.signal import find_peaks
import pickle

mean_spacing_dict = {}
# with open('mean_spacing_cluster.pkl', 'rb') as f:
#     mean_spacing_dict = pickle.load(f)

couplings = [1]

for ratio in np.arange(0.01, 0.21, 0.01):
    temp = np.random.random() + (1/9) / (1 + 1/9)
    mean_spacing = ratio * temp

    energy_grid = np.linspace(0, 10 * temp, 50000)[1:]
    dE = energy_grid[1] - energy_grid[0]

    for mean_x in couplings:
        mean_coupling_strength = np.sqrt(mean_x * mean_spacing / np.pi ** 2)
        total_delay = []

        for i in range(len(total_delay), 1000):
            energy_GOE_mat, W_mat, x_actual = generate_spectrum(num_resonances, mean_spacing, mean_coupling_strength)
            time_delays = wigner_smith_matrix(energy_grid, energy_GOE_mat, np.square(W_mat))
            
            peaks, _ = find_peaks(time_delays.get(), height=10)

            new_energies = np.array([])
            peak_vals = energy_grid[peaks]

            insert = np.array([i * dE / 400 for i in range(400)])

            # if len(peak_vals) > 150:
            #     print("over")
            #     # continue

            if len(peak_vals) > 0:
                for energy in energy_grid:
                    near_peak = False

                    diff = np.abs(energy - peak_vals)
                    if np.min(diff) < 2 * dE:
                        near_peak = True

                    if near_peak:
                        new_energies = np.append(new_energies, energy + insert)
                    else:
                        new_energies = np.append(new_energies, energy)
            else:
                new_energies = energy_grid

            weights = relative_MB_weights(temp, new_energies) * wigner_smith_matrix(new_energies, energy_GOE_mat, np.square(W_mat))
            new_thermal_delay = np.append(new_energies[1:] - new_energies[:-1], dE) * weights

            total_delay.append(np.sum(new_thermal_delay).item())

            if i % 100 == 0:
                with open('mean_spacing_cluster.pkl', 'wb') as f:
                    pickle.dump(mean_spacing_dict, f)

        mean_spacing_dict[mean_spacing.item()] = {'delay': total_delay, 'temp': temp, 'spacing': mean_spacing}

        with open('mean_spacing_cluster.pkl', 'wb') as f:
            pickle.dump(mean_spacing_dict, f)