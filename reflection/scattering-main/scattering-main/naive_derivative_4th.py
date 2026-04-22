import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(suppress=True)

k = 1
dx = 0.01
dx2 = dx ** 2
min_x = -int(1 / k + 1) * 10
max_x = int(1 / k + 1) * 10
wavefunction_e = []
wavefunction_o = []

even = True

def gaussian(mu, sigma):
    return lambda x : -1 / (np.sqrt(2 * np.pi) * sigma) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

def square_well(width, depth):
    return lambda x : -depth if x >= -width and x <= width else 0

V = gaussian(0, 0.01)

def k2(x):
    return k ** 4 - V(x)

def psi(x, index=None):
    global wavefunction
    global even
    if even:
        if x < min_x:
            return np.exp(1j * k * x)
        elif index is not None:
            return wavefunction_e[index]
        else:
            return wavefunction_e[int(round((x - min_x) / dx))]
    else:
        if x < min_x:
            return np.sin(k * x)
        elif index is not None:
            return wavefunction_o[index]
        else:
            return wavefunction_o[int(round((x - min_x) / dx))]


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

    min_x = -1
    max_x = 5

    even = True
    for i in range(0, int((max_x - min_x) / dx) + 1):
        xi = min_x + i * dx
        
        next = 4 * psi(xi - dx, index=i - 1) + (-6 + dx ** 4 * k2(xi - 2 * dx)) * psi(xi - 2 * dx, index=i - 2) + 4 * psi(xi - 3 * dx, index=i - 3) - psi(xi - 4 * dx, index=i - 4)
        print(dx ** 4 * k2(xi - 2 * dx))
        wavefunction_e.append(next)

    # wavefunction = np.array(wavefunction)
    # wavefunction /= max(wavefunction) - min(wavefunction)

    before = [min_x + i * dx for i in range(-int(5 / dx), 0)]
    bpsi2 = [psi(min_x + i * dx) for i in range(-int(5 / dx), 0)]

    after = [min_x + i * dx for i in range(int((max_x - min_x) / dx) + 1)]
    apsi2 = [psi(min_x + i * dx, index=i) for i in range(int((max_x - min_x) / dx) + 1)]

    both = before + after
    pot = [V(x) for x in both]


    x = 0.05
    psi_x = psi(x)
    Dpsi_x = (-psi(x + 2 * dx) + 8 * psi(x + dx) - 8 * psi(x - dx) + psi(x - 2 * dx)) / (12 * dx)

    A = (psi_x * k - 1j * Dpsi_x) * np.exp(-1j * k * x) / (2 * k)
    B = (psi_x * k + 1j * Dpsi_x) * np.exp(1j * k * x) / (2 * k)

    if show:
        print(A, B)
        print(psi_x, Dpsi_x)

    guess = [A * np.exp(1j * k * z) + B * np.exp(-1j * k * z) for z in after]

    if plot:
        plt.plot(after, guess, color='brown', marker='x', markevery=10)
        plt.plot(before, bpsi2, color='red')
        plt.plot(after, apsi2, color='blue')
        plt.plot(both, pot, color='green')
        plt.ylim([min(apsi2) - 1, max(apsi2) + 1])
        plt.show()

    even = False
    for i in range(0, int((max_x - min_x) / dx) + 1):
        xi = min_x + i * dx
        
        next = 4 * psi(xi - dx, index=i - 1) + (-6 + dx ** 4 * k2(xi - 2 * dx)) * psi(xi - 2 * dx, index=i-2) + 4 * psi(xi - 3 * dx, index=i - 3) - psi(xi - 4 * dx, index=i - 4)
        wavefunction_o.append(next)

    # wavefunction = np.array(wavefunction)
    # wavefunction /= max(wavefunction) - min(wavefunction)

    before = [min_x + i * dx for i in range(-int(5 / dx), 0)]
    bpsi2 = [psi(min_x + i * dx) for i in range(-int(5 / dx), 0)]

    after = [min_x + i * dx for i in range(int((max_x - min_x) / dx) + 1)]
    apsi2 = [psi(min_x + i * dx, index=i) for i in range(int((max_x - min_x) / dx) + 1)]

    both = before + after
    pot = [V(x) for x in both]

    x = 1
    psi_x = psi(x)
    Dpsi_x = (-psi(x + 2 * dx) + 8 * psi(x + dx) - 8 * psi(x - dx) + psi(x - 2 * dx)) / (12 * dx)

    C = (psi_x * k - 1j * Dpsi_x) * np.exp(-1j * k * x) / (2 * k)
    D = (psi_x * k + 1j * Dpsi_x) * np.exp(1j * k * x) / (2 * k)

    if show:
        print(C, D)

        print(psi_x, Dpsi_x)

    guess = [C * np.exp(1j * k * z) + D * np.exp(-1j * k * z) for z in after]

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
    return eigenvalues


run(0.01, show=True, plot=True)
exit()

eigenvalues = []
for k_val in np.arange(0.1, 10, 0.1):
    eigen = run(k_val)
    if np.real(eigen[0]) > 0.99:
        eigenvalues.append(eigen[1])
    else:
        eigenvalues.append(eigen[0])
plt.plot(np.arange(0.1, 10, 0.1), np.real(eigenvalues))
plt.show()
plt.plot(np.arange(0.1, 10, 0.1), np.imag(eigenvalues))
plt.show()