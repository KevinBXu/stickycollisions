import numpy as np
import matplotlib.pyplot as plt

k = 0.01
dx = 0.01
dx2 = dx ** 2
min_x = -int(1 / k + 1) * 10
max_x = int(1 / k + 1) * 10
wavefunction_e = []
wavefunction_o = []

def V(x):
    # return 0
    # if x < 10 and x > -10:
    #     return -10
    # else:
    #     return 0
    mu = 0
    sigma = 0.01
    return -1 / (np.sqrt(2 * np.pi) * sigma) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

def k2(x):
    return k ** 2 - V(x)

def psi(x, index=None):
    global wavefunction
    global even
    if even:
        if x > max_x:
            return np.cos(k * x)
        elif index is not None:
            return wavefunction_e[index]
        else:
            return wavefunction_e[int(round((x - min_x) / dx))]
    else:
        if x > max_x:
            return np.sin(k * x)
        elif index is not None:
            return wavefunction_o[index]
        else:
            return wavefunction_o[int(round((x - min_x) / dx))]

even = True
for i in range(int((max_x - min_x) / dx), -1, -1):
    xi = min_x + i * dx
    
    a = psi(xi - dx, index=i - 1) * (2 + 10 / 12 * dx2 * k2(xi - dx))
    b = psi(xi - 2 * dx, index=i - 2) * (1 - 1 / 12 * dx2 * k2(xi - 2 * dx))
    c = 1 - 1 / 12 * dx2 * k2(xi)
    next = (a - b) / c
    wavefunction_e.append(next)

# wavefunction = np.array(wavefunction)
# wavefunction /= max(wavefunction) - min(wavefunction)

before = [min_x + i * dx for i in range(-int(5 / dx), 0)]
bpsi2 = [psi(min_x + i * dx) for i in range(-int(5 / dx), 0)]

after = [min_x + i * dx for i in range(int((max_x - min_x) / dx) + 1)]
apsi2 = [psi(min_x + i * dx, index=i) for i in range(int((max_x - min_x) / dx) + 1)]

both = before + after
pot = [V(x) for x in both]


x = 5
psi_x = psi(x)
Dpsi_x = (-psi(x + 2 * dx) + 8 * psi(x + dx) - 8 * psi(x - dx) + psi(x - 2 * dx)) / (12 * dx)

A = (psi_x * k - 1j * Dpsi_x) * np.exp(-1j * k * x) / (2 * k)
B = (psi_x * k + 1j * Dpsi_x) * np.exp(1j * k * x) / (2 * k)

print(A, B)

print(psi_x, Dpsi_x)

guess = [A * np.exp(1j * k * z) + B * np.exp(-1j * k * z) for z in after]

plt.plot(after, guess, color='brown', marker='x', markevery=1000)
plt.plot(before, bpsi2, color='red')
plt.plot(after, apsi2, color='blue')
plt.plot(both, pot, color='green')
plt.ylim([min(apsi2) - 1, max(apsi2) + 1])
plt.show()