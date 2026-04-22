import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import lu

np.set_printoptions(suppress=True)
np.set_printoptions(linewidth=np.inf)

dim = 4
cosines = []
sines = []

dx = 0.001
dx2 = dx ** 2
wavefunction = []
func = []

even = True

def gaussian(mu, sigma):
    return lambda x : -1 / (np.sqrt(2 * np.pi) * sigma) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

def square_well(width, depth):
    return lambda x : -depth if x >= -width and x <= width else 0

def zero():
    return lambda x : 0

V = gaussian(0, 0.5)

def psi(x, index=None):
    global wavefunction
    global func

    if x < min_x:
        return np.array(func(x))
    elif index is not None:
        return wavefunction[index]
    else:
        return wavefunction[int(round((x - min_x) / dx))]
        
def f(xi, psi_i):
    vec = np.array([1j for _ in range(dim)])
    for i in range(dim - 1):
        vec[i] = psi_i[i + 1]
    vec[dim - 1] = ((-1j) ** dim) * (k ** dim - V(xi)) * psi_i[0]
    return vec


def run(k_val, show=False, plot=False):
    global wavefunction
    global func
    global min_x
    global max_x
    global k
    global sines
    global cosines
    k = k_val
    roots = [1j * k, -1j * k, k, -k]

    min_x = -5
    max_x = 15

    wavefunction = []

    coef = [1, 1, (-0.40000048951115086-0j), 0]

    func = lambda x : [sum([coef[j] * roots[j] ** i * np.exp(roots[j] * x) for j in range(dim)]) for i in range(dim)]

    for i in range(-1, int((max_x - min_x) / dx) + 1):
        xi = min_x + i * dx
        psi_i = psi(xi, index=i)
        
        k1 = dx * f(xi, psi_i)
        k2 = dx * f(xi + dx, psi_i + k1)

        next = psi_i + (k1 + k2) / 2

        wavefunction.append(next)

    before = [min_x + i * dx for i in range(-int(10 / dx), 0)]
    bpsi2 = [np.abs(psi(xi)[0]) ** 2 for xi in before]

    after = [min_x + i * dx for i in range(int((max_x - min_x) / dx) + 1)]
    apsi2 = [np.abs(psi(xi, index=i)[0]) ** 2 for i in range(int((max_x - min_x) / dx) + 1)]

    both = before + after
    pot = [V(x) for x in both]

    x = 0.5


    a = np.full((dim, dim), 1j)
    for der in range(dim):
        for (i, mom) in enumerate(roots):
            a[der][i] = mom ** der * np.exp(mom * x)
    b = psi(x)
    print("Matrices")
    print(a)
    print(b)
    x = np.linalg.solve(a, b)
    print("Solution")
    print(x)




    # guess = [np.abs(A * np.exp(1j * k * z) + B * np.exp(-1j * k * z)) ** 2 for z in after]

    if plot:
        # plt.plot(after, guess, color='brown', marker='x', markevery=10)
        plt.plot(before, bpsi2, color='red')
        plt.plot(after, apsi2, color='blue')
        plt.plot(both, pot, color='green')
        plt.title("Quartic")
        plt.xlabel(r"$x$")
        plt.ylabel(r"$\psi(x)$")
        plt.ylim([-10, 10])
        plt.show()


run(1, show=True, plot=True)