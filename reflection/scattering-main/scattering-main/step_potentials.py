import numpy as np
import matplotlib.pyplot as plt

plot = False
Ls = [2, 1, 0.1, 0.001]
Es = np.arange(0.0001, 0.005, 0.0001)
order = 1/5

def gaussian(sigma):
    # return 0
    # if x < 10 and x > -10:
    #     return -10
    # else:
    #     return 0
    mu = 0
    return lambda x : -1 / (np.sqrt(2 * np.pi) * sigma) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

def step(x):
    if x < 0:
        return -1
    elif x < L: 
        return -1 + (x / L) ** order 
    else:
        return 0

V = lambda x : step(x)

for order in [0.2]:
    for L in Ls:
        r = []
        for E in Es:
            E = E 
            k = np.sqrt(E)
            k_prime = np.sqrt(E - V(-1))
        
            dx = 0.01
            dx2 = dx ** 2
            min_x = -20
            max_x = 20
            wavefunction_e = []
            wavefunction_o = []

            def k2(x):
                return -E + V(x)

            def psi(x, index=None):
                global wavefunction
                global even
                if even:
                    if x < min_x:
                        return np.cos(k_prime * x)
                    elif index is not None:
                        return wavefunction_e[index]
                    else:
                        return wavefunction_e[int(round((x - min_x) / dx))]
                else:
                    if x < min_x:
                        return np.sin(k_prime * x)
                    elif index is not None:
                        return wavefunction_o[index]
                    else:
                        return wavefunction_o[int(round((x - min_x) / dx))]

            even = True
            for i in range(0, int((max_x - min_x) / dx) + 1):
                xi = min_x + i * dx
                
                a = psi(xi - dx, index=i - 1) * (2 + 10 / 12 * dx2 * k2(xi - dx))
                b = psi(xi - 2 * dx, index=i - 2) * (1 - 1 / 12 * dx2 * k2(xi - 2 * dx))
                c = 1 - 1 / 12 * dx2 * k2(xi)
                next = (a - b) / c
                wavefunction_e.append(next)

            # wavefunction = np.array(wavefunction)
            # wavefunction /= max(wavefunction) - min(wavefunction)

            before = [min_x + i * dx for i in range(-int(5 / dx), int((max_x - min_x) / dx) + 1)]
            bpsi2 = [np.cos(k_prime * x) for x in before]

            after = [min_x + i * dx for i in range(int((max_x - min_x) / dx) + 1)]
            apsi2 = [psi(min_x + i * dx, index=i) for i in range(int((max_x - min_x) / dx) + 1)]

            both = before
            pot = [V(x) for x in both]


            x = 5
            psi_x = psi(x)
            Dpsi_x = (-psi(x + 2 * dx) + 8 * psi(x + dx) - 8 * psi(x - dx) + psi(x - 2 * dx)) / (12 * dx)

            A = (psi_x * k - 1j * Dpsi_x) * np.exp(-1j * k * x) / (2 * k)
            B = (psi_x * k + 1j * Dpsi_x) * np.exp(1j * k * x) / (2 * k)

            # print(A, B)

            # print(psi_x, Dpsi_x)

            guess = [A * np.exp(1j * k * z) + B * np.exp(-1j * k * z) for z in after]

            if plot:
                plt.title("Even" if even else "Odd")
                plt.plot(after, guess, color='brown', marker='x', markevery=10)
                plt.plot(before, bpsi2, color='red')
                plt.plot(after, apsi2, color='blue')
                plt.plot(both, pot, color='green')
                plt.ylim([min(apsi2) - 1, max(apsi2) + 1])
                plt.show()

            even = False
            for i in range(0, int((max_x - min_x) / dx) + 1):
                xi = min_x + i * dx
                
                a = psi(xi - dx, index=i - 1) * (2 + 10 / 12 * dx2 * k2(xi - dx))
                b = psi(xi - 2 * dx, index=i - 2) * (1 - 1 / 12 * dx2 * k2(xi - 2 * dx))
                c = 1 - 1 / 12 * dx2 * k2(xi)
                next = (a - b) / c
                wavefunction_o.append(next)

            # wavefunction = np.array(wavefunction)
            # wavefunction /= max(wavefunction) - min(wavefunction)

            before = [min_x + i * dx for i in range(-int(5 / dx), int((max_x - min_x) / dx) + 1)]
            bpsi2 = [np.sin(k_prime * x) for x in before]

            after = [min_x + i * dx for i in range(int((max_x - min_x) / dx) + 1)]
            apsi2 = [psi(min_x + i * dx, index=i) for i in range(int((max_x - min_x) / dx) + 1)]

            both = before
            pot = [V(x) for x in both]

            x = 5
            psi_x = psi(x)
            Dpsi_x = (-psi(x + 2 * dx) + 8 * psi(x + dx) - 8 * psi(x - dx) + psi(x - 2 * dx)) / (12 * dx)

            C = (psi_x * k - 1j * Dpsi_x) * np.exp(-1j * k * x) / (2 * k)
            D = (psi_x * k + 1j * Dpsi_x) * np.exp(1j * k * x) / (2 * k)

            # print(C, D)

            # print(psi_x, Dpsi_x)

            guess = [C * np.exp(1j * k * z) + D * np.exp(-1j * k * z) for z in after]

            if plot:
                plt.title("Even" if even else "Odd")
                plt.plot(after, guess, color='brown', marker='x', markevery=10)
                plt.plot(before, bpsi2, color='red')
                plt.plot(after, apsi2, color='blue')
                plt.plot(both, pot, color='green')
                plt.ylim([min(apsi2) - 1, max(apsi2) + 1])
                plt.show()

            left = np.array([[A, C], [1 / 2, 1j / 2]])
            right = np.array([[1 / 2, -1j / 2], [B, D]])
            S = left @ np.linalg.inv(right)
            # print(S)

            eigenvalues, eigenvectors = np.linalg.eig(S)
            # print("Reflection", np.linalg.norm((S @ [1, 0])[1]))
            r.append(np.linalg.norm((S @ [1, 0])[1]))

        plt.plot(Es, r, label="dV/dx = " + str(L) + "d=" + str(order))

plt.title("Reflection Coefficient for d = " + str(order))
plt.xlabel("E - V")
plt.ylabel(f"R^2")
plt.legend()
plt.show()