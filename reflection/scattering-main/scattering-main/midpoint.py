import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(suppress=True)

cosines = [np.cos, lambda x : -np.sin(x), lambda x : -np.cos(x), np.sin]
sines = [np.sin, np.cos, lambda x : -np.sin(x), lambda x : -np.cos(x)]

dim = 4

dx = 0.01
dx2 = dx ** 2
wavefunction_e = []
wavefunction_o = []

even = True

def gaussian(mu, sigma):
    return lambda x : -1 / (np.sqrt(2 * np.pi) * sigma) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

def square_well(width, depth):
    return lambda x : -depth if x >= -width and x <= width else 0

def zero():
    return lambda x : 0

V = square_well(0.01, 0.01)

def psi(x, index=None):
    global wavefunction
    global even
    
    if even:
        if x < min_x:
            return np.array([k ** i * cosines[i](k * x) for i in range(dim)])
        elif index is not None:
            return wavefunction_e[index]
        else:
            return wavefunction_e[int(round((x - min_x) / dx))]
    else:
        if x < min_x:
            return np.array([k ** i * sines[i](k * x) for i in range(dim)])
        elif index is not None:
            return wavefunction_o[index]
        else:
            return wavefunction_o[int(round((x - min_x) / dx))]
        
def f(xi, psi_i):
    vec = np.array([1j for _ in range(dim)])
    for i in range(dim - 1):
        vec[i] = psi_i[i + 1]
    vec[dim - 1] = ((-1j) ** dim) * (k ** 2 - V(xi)) * psi_i[0]
    return vec


def run(k_val, show=False, plot=False):
    global even
    global wavefunction_e
    global wavefunction_o
    global min_x
    global max_x
    global k
    wavefunction_e = []
    wavefunction_o = []
    k = k_val

    min_x = -10
    max_x = 10

    even = True
    for i in range(-1, int((max_x - min_x) / dx) + 1):
        xi = min_x + i * dx
        psi_i = psi(xi, index=i)
        
        k1 = dx * f(xi, psi_i)
        k2 = dx * f(xi + dx, psi_i + k1)

        next = psi_i + (k1 + k2) / 2

        wavefunction_e.append(next)

    # wavefunction = np.array(wavefunction)
    # wavefunction /= max(wavefunction) - min(wavefunction)

    before = [min_x + i * dx for i in range(-int(5 / dx), 0)]
    bpsi2 = [np.abs(psi(min_x + i * dx)[0]) ** 2 for i in range(-int(5 / dx), 0)]

    after = [min_x + i * dx for i in range(int((max_x - min_x) / dx) + 1)]
    apsi2 = [np.abs(psi(min_x + i * dx, index=i)[0]) ** 2 for i in range(int((max_x - min_x) / dx) + 1)]

    both = before + after
    pot = [V(x) for x in both]


    x = 1
    psi_x = psi(x)[0]
    Dpsi_x = (-psi(x + 2 * dx) + 8 * psi(x + dx) - 8 * psi(x - dx) + psi(x - 2 * dx))[0] / (12 * dx)

    A = (psi_x * k - 1j * Dpsi_x) * np.exp(-1j * k * x) / (2 * k)
    B = (psi_x * k + 1j * Dpsi_x) * np.exp(1j * k * x) / (2 * k)

    if show:
        print(A, B)
        print(psi_x, Dpsi_x)

    guess = [np.abs(A * np.exp(1j * k * z) + B * np.exp(-1j * k * z)) ** 2 for z in after]

    if plot:
        plt.plot(after, guess, color='brown', marker='x', markevery=10)
        plt.plot(before, bpsi2, color='red')
        plt.plot(after, apsi2, color='blue')
        plt.plot(both, pot, color='green')
        plt.ylim([min(apsi2) - 1, max(apsi2) + 1])
        plt.show()

    even = False
    for i in range(-1, int((max_x - min_x) / dx) + 1):
        xi = min_x + i * dx
        psi_i = psi(xi, index=i)
        
        k1 = dx * f(xi, psi_i)
        k2 = dx * f(xi + dx, psi_i + k1)

        next = psi_i + (k1 + k2) / 2
        wavefunction_o.append(next)

    # wavefunction = np.array(wavefunction)
    # wavefunction /= max(wavefunction) - min(wavefunction)

    before = [min_x + i * dx for i in range(-int(5 / dx), 0)]
    bpsi2 = [np.abs(psi(min_x + i * dx)[0]) ** 2 for i in range(-int(5 / dx), 0)]

    after = [min_x + i * dx for i in range(int((max_x - min_x) / dx) + 1)]
    apsi2 = [np.abs(psi(min_x + i * dx, index=i)[0]) ** 2 for i in range(int((max_x - min_x) / dx) + 1)]

    both = before + after
    pot = [V(x) for x in both]

    x = 1
    psi_x = psi(x)[0]
    Dpsi_x = (-psi(x + 2 * dx) + 8 * psi(x + dx) - 8 * psi(x - dx) + psi(x - 2 * dx))[0] / (12 * dx)

    C = (psi_x * k - 1j * Dpsi_x) * np.exp(-1j * k * x) / (2 * k)
    D = (psi_x * k + 1j * Dpsi_x) * np.exp(1j * k * x) / (2 * k)

    if show:
        print(C, D)

        print(psi_x, Dpsi_x)

    guess = [np.abs(C * np.exp(1j * k * z) + D * np.exp(-1j * k * z)) ** 2 for z in after]

    if plot:
        plt.plot(after, guess, color='brown', marker='x', markevery=10)
        plt.plot(before, bpsi2, color='red')
        plt.plot(after, apsi2, color='blue')
        plt.plot(both, pot, color='green')
        plt.ylim([min(apsi2) - 1, max(apsi2) + 1])
        plt.show()

    left = np.array([[A, C], [1 / 2, 1j / 2]])
    right = np.array([[1 / 2, -1j / 2], [B, D]])
    S = left @ np.linalg.inv(right)
    inverse = np.array([[1, 1], [1, -1]])


    eigenvalues, eigenvectors = np.linalg.eig(S)

    if show:
        print('Scattering Matrix', S)
        print('Even Odd Basis', inverse @ S @ inverse / 2)
        print('eigenvalues', eigenvalues, np.abs(eigenvalues[0]), np.abs(eigenvalues[1]))
        print("det", np.abs(np.linalg.det(S)))
        print('size', len(wavefunction_e))
    return eigenvalues


run(1, show=True, plot=True)