import cupy as np
import math
import matplotlib.pyplot as plt
import scipy as sc
import numpy

# parameters
e_mass = 5.48579909-4     # electron mass in amu
t0 = 315775               # a.u. to K
tau0 = 2.418884326509e-17 # a.u. to sec
 
# universal truths, in SI units
kB = 1.3806e-23
hbar = 1.05457e-34
AMU = 1.6605e-27
Hartree = 4.359744e-18
auLength = 5.29177e-11

# collision parameters
abar = (np.pi * 2 ** (-3/2) / math.gamma(5/4) / math.gamma(1/2)) ** 2

def generate_spectrum(num_resonances, mean_spacing, mean_coupling_strength, num_open_channels = 1):
    # assert(num_resonances % 2 == 1)
    # midpoint = num_resonances // 2
    # energy_max = num_resonances / 2 * mean_spacing
    # energy_min = -energy_max
    
    # rng
    rng = numpy.random.default_rng()

    energy_GOE_mat = np.zeros(num_resonances * len(mean_coupling_strength))
    W_mat = np.zeros(num_resonances * len(mean_coupling_strength))

    for j in range(len(mean_coupling_strength)):
        uniform_dist = np.array(rng.uniform(0, 1, num_resonances - 1))

        # transform uniform into a Wigner-Dyson distribution
        nearest_neighbor = mean_spacing * np.sqrt(-4 / np.pi * np.log(1 - uniform_dist))

        starting_index = j * num_resonances
        energy_GOE_mat[starting_index] = -np.sum(nearest_neighbor) / 2
        for i in range(num_resonances - 1):
            energy_GOE_mat[starting_index + i + 1] = energy_GOE_mat[starting_index + i] + nearest_neighbor[i]
        
        for i in range(num_resonances):
            if energy_GOE_mat[i] == 0:
                print("mayday")

        # generate coupling matrix
        for mu in range(num_resonances):
            W_mat[starting_index + mu] = np.array(rng.normal(0, mean_coupling_strength[j].get()))

    # find effective hamiltonian
    # Heff = energy_GOE_mat - 1j * np.pi * (W_mat @ W_mat.T) 

    # actual coupling constant
    x_actual = np.pi ** 2 / mean_spacing * np.var(W_mat)

    return energy_GOE_mat, W_mat, x_actual

def generate_spectrum_sectors_funnel(num_resonances, mean_spacing, mean_coupling_strength, num_open_channels = 1):
    # assert(num_resonances % 2 == 1)
    # midpoint = num_resonances // 2
    # energy_max = num_resonances / 2 * mean_spacing
    # energy_min = -energy_max
    
    # rng
    rng = numpy.random.default_rng()

    energy_GOE_mat = np.zeros((num_resonances * len(mean_coupling_strength), num_resonances * len(mean_coupling_strength)))
    W_mat = np.zeros(num_resonances * len(mean_coupling_strength))

    for j in range(len(mean_coupling_strength)):
        uniform_dist = np.array(rng.uniform(0, 1, num_resonances - 1))

        # transform uniform into a Wigner-Dyson distribution
        nearest_neighbor = mean_spacing * np.sqrt(-4 / np.pi * np.log(1 - uniform_dist))

        starting_index = j * num_resonances
        energy_GOE_mat[starting_index, starting_index] = -np.sum(nearest_neighbor) / 2
        for i in range(num_resonances - 1):
            energy_GOE_mat[starting_index + i + 1, starting_index + i + 1] = energy_GOE_mat[starting_index + i, starting_index + i] + nearest_neighbor[i]
        
        for i in range(num_resonances):
            if energy_GOE_mat[i, i] == 0:
                print("mayday")

        # generate coupling matrix
        if j == 0:
            for mu in range(num_resonances):
                W_mat[starting_index + mu] = np.array(rng.normal(0, mean_coupling_strength[j].get()))

    for i in range(1, len(mean_coupling_strength)):
        starting_index = i * num_resonances
        couplings = np.abs(np.array(rng.normal(0, mean_coupling_strength[i].get(), size=(num_resonances, num_resonances))))
        energy_GOE_mat[0:num_resonances, starting_index:starting_index + num_resonances] = couplings
        energy_GOE_mat[starting_index:starting_index + num_resonances, 0:num_resonances] = couplings.T


    # print(energy_GOE_mat)
    # print(W_mat)
    if len(mean_coupling_strength) > 1 and mean_coupling_strength[1] != 0:
        eigenvalues, eigenvectors = np.linalg.eigh(energy_GOE_mat) 
        energy_GOE_mat = np.conj(eigenvectors.T) @ energy_GOE_mat @ eigenvectors
        W_mat = np.conj(eigenvectors.T) @ W_mat

    energy_GOE_mat = np.diag(energy_GOE_mat)

    # find effective hamiltonian
    # Heff = energy_GOE_mat - 1j * np.pi * (W_mat @ W_mat.T) 

    # actual coupling constant
    x_actual = np.pi ** 2 / mean_spacing * np.var(W_mat)

    return energy_GOE_mat, W_mat, x_actual

def generate_spectrum_sectors_cascade(num_resonances, mean_spacing, mean_coupling_strength, num_open_channels = 1):
    # assert(num_resonances % 2 == 1)
    # midpoint = num_resonances // 2
    # energy_max = num_resonances / 2 * mean_spacing
    # energy_min = -energy_max
    
    # rng
    rng = numpy.random.default_rng()

    energy_GOE_mat = np.zeros((num_resonances * len(mean_coupling_strength), num_resonances * len(mean_coupling_strength)))
    W_mat = np.zeros(num_resonances * len(mean_coupling_strength))

    for j in range(len(mean_coupling_strength)):
        uniform_dist = np.array(rng.uniform(0, 1, num_resonances - 1))

        # transform uniform into a Wigner-Dyson distribution
        nearest_neighbor = mean_spacing * np.sqrt(-4 / np.pi * np.log(1 - uniform_dist))

        starting_index = j * num_resonances
        energy_GOE_mat[starting_index, starting_index] = -np.sum(nearest_neighbor) / 2
        for i in range(num_resonances - 1):
            energy_GOE_mat[starting_index + i + 1, starting_index + i + 1] = energy_GOE_mat[starting_index + i, starting_index + i] + nearest_neighbor[i]
        
        for i in range(num_resonances):
            if energy_GOE_mat[i, i] == 0:
                print("mayday")

        # generate coupling matrix
        if j == 0:
            for mu in range(num_resonances):
                W_mat[starting_index + mu] = np.array(rng.normal(0, mean_coupling_strength[j].get()))

    for i in range(0, len(mean_coupling_strength) - 1):
        index1 = i * num_resonances
        index2 = (i+1) * num_resonances
        couplings = np.abs(np.array(rng.normal(0, mean_coupling_strength[i].get(), size=(num_resonances, num_resonances))))
        energy_GOE_mat[index2:index2 + num_resonances, index1:index1 + num_resonances] = couplings
        energy_GOE_mat[index1:index1 + num_resonances, index2:index2 + num_resonances] = couplings.T


    # print(energy_GOE_mat)
    # print(W_mat)
    if len(mean_coupling_strength) > 1 and mean_coupling_strength[1] != 0:
        eigenvalues, eigenvectors = np.linalg.eigh(energy_GOE_mat) 
        energy_GOE_mat = np.conj(eigenvectors.T) @ energy_GOE_mat @ eigenvectors
        W_mat = np.conj(eigenvectors.T) @ W_mat

    energy_GOE_mat = np.diag(energy_GOE_mat)

    # find effective hamiltonian
    # Heff = energy_GOE_mat - 1j * np.pi * (W_mat @ W_mat.T) 

    # actual coupling constant
    x_actual = np.pi ** 2 / mean_spacing * np.var(W_mat)

    return energy_GOE_mat, W_mat, x_actual

def wigner_smith_matrix(collision_energy, H_goe, W_squared):
    # values for calculation

    k = np.sqrt(collision_energy)
    A_mqdt = abar * k
    A_prime = abar / 2 / k
    G_mqdt = (1 / 3 - abar ** 2) * np.square(k)
    G_prime = (1 / 3 - abar ** 2)

    collision_energy = collision_energy.reshape((-1, 1))
    energy_diff = collision_energy - H_goe
    ratio = W_squared / energy_diff

    Y_mat = -np.pi * np.sum(ratio, axis=1)
    Y_prime = np.pi * np.sum(ratio / energy_diff, axis=1)

    K_mat = A_mqdt * Y_mat / (1 + G_mqdt * Y_mat)
    K_prime = (A_prime * Y_mat + A_mqdt * Y_prime) / (1 + G_mqdt * Y_mat) - A_mqdt * Y_mat * (G_prime * Y_mat + G_mqdt * Y_prime) / (1 + G_mqdt * Y_mat) ** 2

    Q_mat = 2 * (K_prime / (1 + K_mat ** 2))
    return Q_mat

def S_matrix_alt(collision_energy, H_goe, W_mat):
    num_res = len(H_goe)

    D = np.diag(collision_energy - H_goe) + 1j * np.pi * np.square(W_mat)

    S = 1 - 2j * np.pi * np.inner(W_mat, numpy.linalg.inv(D) @ W_mat)

    return S


def wigner_smith_matrix_single(collision_energy, H_goe, W_squared):
    # values for calculation
    k = np.sqrt(collision_energy)
    A_mqdt = abar * k
    A_prime = abar / 2 / k
    G_mqdt = (1 / 3 - abar ** 2) * np.square(k)
    G_prime = (1 / 3 - abar ** 2)
    eta_mqdt = -abar * k
    eta_prime = -abar / 2 / k 

    energy_diff = collision_energy - H_goe
    ratio = W_squared / energy_diff

    Y_mat = -np.pi * np.sum(ratio)
    Y_prime = np.pi * np.sum(ratio / energy_diff)

    K_mat = A_mqdt * Y_mat / (1 + G_mqdt * Y_mat)
    K_prime = (A_prime * Y_mat + A_mqdt * Y_prime) / (1 + G_mqdt * Y_mat) - A_mqdt * Y_mat * (G_prime * Y_mat + G_mqdt * Y_prime) / (1 + G_mqdt * Y_mat) ** 2

    Q_mat = 2 * (K_prime / (1 + K_mat ** 2))
    
    return Q_mat

def S_matrix(collision_energy, H_goe, W_squared):
    # values for calculation
    k = np.sqrt(collision_energy).reshape(-1)
    A_mqdt = abar * k
    G_mqdt = (1 / 3 - abar ** 2) * np.square(k)
    eta_mqdt = -abar * k

    collision_energy = np.reshape(collision_energy, (len(collision_energy), 1))
    energy_diff = collision_energy - H_goe
    ratio = W_squared / energy_diff

    Y_mat = -np.pi * np.sum(ratio, axis=1)

    K_mat = A_mqdt * Y_mat / (1 + G_mqdt * Y_mat)
    S_mat = (1 + 1j * K_mat) / (1 - 1j * K_mat)

    return S_mat    


def relative_MB_weights(temp, energy):
    weight = 2 * np.sqrt(energy / np.pi) * np.exp(-energy / temp) / temp ** (3/2) 
    return weight

def thermal_time_delay(temp, energy_GOE_mat, W_mat, num_energies):
    energy_grid = np.linspace(0, 10 * temp, num=num_energies)[1:]
    grid_spacing = energy_grid[1] - energy_grid[0]

    relative_time_delay = grid_spacing * relative_MB_weights(temp, energy_grid) * wigner_smith_matrix(energy_grid, energy_GOE_mat, np.square(W_mat))

    return np.sum(relative_time_delay)
    