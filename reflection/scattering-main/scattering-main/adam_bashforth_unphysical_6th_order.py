import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import lu
from scipy import special
import pickle
import sys

orig_stdout = sys.stdout
f = open('log.log', 'w')
sys.stdout = f

np.set_printoptions(suppress=True)
np.set_printoptions(linewidth=np.inf)
np.set_printoptions(precision=3)

dim = 6
cosines = []
sines = []
coefficients = []

dx = 0.001
dx2 = dx ** 2
wavefunction = []
func = []
c = None
h = None

even = True

def gaussian(sigma):
    return lambda x : -1 / (np.sqrt(2 * np.pi) * sigma) * np.exp(-0.5 * (x / sigma) ** 2)

def gaussian_1(sigma):
    return lambda x : (x / sigma ** 2) / (np.sqrt(2 * np.pi) * sigma) * np.exp(-0.5 * (x / sigma) ** 2)

def der_gaussian(sigma, n):
    c = (-1) ** (n + 1) / (np.sqrt(2 * np.pi) * sigma) * (1 / 2 / sigma ** 2) ** (n / 2)
    h = special.hermite(n)
    return lambda x : c * np.exp(-0.5 * (x / sigma) ** 2) * h(x / sigma / np.sqrt(2))

def square_well(width, depth):
    return lambda x : -depth if x >= -width and x <= width else 0

def lorentzian(width):
    return lambda x : -width / (x ** 2 + width ** 2)

def zero():
    return lambda x : 0

def f():
    return lambda x : -np.exp(-x ** 2 / 4) * (x ** 2 - 2)

V = gaussian(0.2)

def psi(x, index=None):
    global wavefunction
    global func

    if x < min_x:
        return np.array(func(x), dtype=np.cdouble)
    elif index is not None:
        return wavefunction[index]
    else:
        return wavefunction[int(round((x - min_x) / dx))]
        
def f(xi, psi_i):
    vec = np.array([1j for _ in range(dim)], dtype=np.cdouble)
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
    roots = np.array([1j * k * np.exp(2 * np.pi * 1j * l / dim) for l in range(dim)], dtype=np.cdouble)
    print(roots)

    min_x = -4
    max_x = 5
    x = 5

    for momentum in roots:
        wavefunction = []

        func = lambda x : np.array([momentum ** i * np.exp(momentum * x) for i in range(dim)], dtype=np.cdouble)

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

        a = np.full((dim, dim), 1j, dtype=np.cdouble)
        for der in range(dim):
            for (i, mom) in enumerate(roots):
                a[der][i] = mom ** der * np.exp(mom * x)
        b = psi(x)
        # print("Matrices")
        # print(a)
        # print(b)
        sol = np.linalg.solve(a, b)
        print("Solution")
        print(sol)

        coefficients.append(sol)


        # guess = [np.abs(A * np.exp(1j * k * z) + B * np.exp(-1j * k * z)) ** 2 for z in after]

        if plot:
            plot = False
            # plt.plot(after, guess, color='brown', marker='x', markevery=10)
            plt.plot(before, bpsi2, color='red')
            plt.plot(after, apsi2, color='blue')
            plt.plot(both, pot, color='green')
            plt.ylim([-1, 10])
            plt.show()

    
    S = np.array(coefficients, dtype=np.cdouble).T
    with open('S-matrix.pkl', 'wb') as file:
        pickle.dump(S, file)

    def compute(v, even):
        vector = 1j * np.zeros(dim, dtype=np.cdouble)
        vector[0] = 1
        if even:
            vector[dim // 2] = 1
        else:
            vector[dim // 2] = -1
        for i in range(1, dim // 2):
            vector[i] = v[i - 1]
        return S @ vector

    coeffs = S[1:dim // 2, 1:dim // 2]

    even_output = -compute(np.zeros(dim // 2, dtype=np.cdouble), True)[1:dim // 2]
    even_sol = np.linalg.solve(coeffs, even_output)
    x1 = compute(even_sol, True)

    odd_output = -compute(np.zeros(dim // 2, dtype=np.cdouble), False)[1:dim // 2]
    coeffs = S[1:dim // 2, 1:dim // 2]
    odd_sol = np.linalg.solve(coeffs, odd_output)
    x2 = compute(odd_sol, False)

    print("x1", x1)
    print("x2", x2)

    left = np.array([[x1[0], x2[0]], [1, -1]], dtype=np.cdouble)
    right = np.array([[1, 1], [x1[dim // 2], x2[dim // 2]]], dtype=np.cdouble)
    S2 = left @ np.linalg.inv(right)
    # print(S2 @ np.array([1, x1[2]]))
    print("S2")
    print(S2)
    print(S2 @ [1, 1])
    eigenvalues, eigenvectors = np.linalg.eig(S2)
    eigenvalues = np.conj(eigenvalues)
    print("Eigenvalues")
    print(eigenvalues)
    print("Eigenvectors")
    print(eigenvectors)


run(0.01, show=True, plot=True)