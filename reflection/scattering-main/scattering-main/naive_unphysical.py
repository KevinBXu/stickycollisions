import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import lu
import math

np.set_printoptions(suppress=True)
np.set_printoptions(linewidth=np.inf)

dim = 2
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

V = zero()

def psi(x, index=None):
    global wavefunction
    global func

    if x < min_x:
        return np.array([np.exp(1j * k * x) for i in range(1)])
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

def k2(x):
    return k ** dim - V(x)

def run(k_val, show=False, plot=False):
    global wavefunction
    global func
    global min_x
    global max_x
    global k
    global sines
    global cosines
    k = k_val

    min_x = -5
    max_x = 10

    mid = dim // 2
    for momentum in [1j * k]: #, -1j * k, k, -k
        wavefunction = []

        func = [lambda x : (momentum ** i) * np.exp(momentum * x) for i in range(dim)]

        for i in range(0, int((max_x - min_x) / dx) + 1):
            xi = min_x + i * dx
            
            next = (-1j) ** dim * dx ** dim * k2(xi - mid * dx) * psi(xi - mid * dx, index=i - mid)
            for j in range(0, dim):
                next -= (-1) ** (j + dim) * math.comb(dim, j) * psi(xi - (dim - j) * dx, index=i - dim + j)

            wavefunction.append(next)

        before = [min_x + i * dx for i in range(-int(10 / dx), 0)]
        bpsi2 = [np.abs(psi(xi)[0]) ** 2 for xi in before]

        after = [min_x + i * dx for i in range(int((max_x - min_x) / dx) + 1)]
        apsi2 = [np.abs(psi(xi, index=i)[0]) ** 2 for i in range(int((max_x - min_x) / dx) + 1)]

        both = before + after
        pot = [V(x) for x in both]

        x = 0.5


        # a = np.full((dim, dim), 1j)
        # for der in range(dim):
        #     for (i, mom) in enumerate([1j * k, -1j * k, k, -k]):
        #         a[der][i] = mom ** der * np.exp(mom * x)
        # b = psi(x)
        # print("Matrices")
        # print(a)
        # print(b)
        # x = np.linalg.solve(a, b)
        # print("Solution")
        # print(x)




        # guess = [np.abs(A * np.exp(1j * k * z) + B * np.exp(-1j * k * z)) ** 2 for z in after]

        if plot:
            # plt.plot(after, guess, color='brown', marker='x', markevery=10)
            plt.plot(before, bpsi2, color='red')
            plt.plot(after, apsi2, color='blue')
            plt.plot(both, pot, color='green')
            plt.ylim([-1, 10])
            plt.show()

    return 


run(1, show=True, plot=True)