import matplotlib.pyplot as plt
import math
import numpy as np
from scipy import special

def der_gaussian(sigma, n):
    c = (-1) ** (n + 1) / (np.sqrt(2 * np.pi) * sigma) * (1 / 2 / sigma ** 2) ** (n / 2)
    h = special.hermite(n)
    return lambda x : c * np.exp(-0.5 * (x / sigma) ** 2) * h(x / sigma / np.sqrt(2))

def gaussian(sigma):
    return lambda x : -1 / (np.sqrt(2 * np.pi) * sigma) * np.exp(-0.5 * (x / sigma) ** 2)

def gaussian_1(sigma):
    return lambda x : -(x / (sigma ** 2)) * gaussian(sigma)(x)

domain = np.arange(-10, 10, 0.01)

plt.plot(domain, [gaussian(0.005)(x) for x in domain])
plt.plot(domain, [der_gaussian(0.005, 0)(x) for x in domain])
plt.show()

plt.plot(domain, [gaussian_1(0.005)(x) for x in domain])
plt.plot(domain, [der_gaussian(0.005, 1)(x) for x in domain])
plt.show()
plt.plot(domain, np.array([gaussian_1(0.7)(x) for x in domain]) / np.array([der_gaussian(0.7, 1)(x) for x in domain]))
plt.show()

exit()
for i in range(4):
	g = der_gaussian(0.1, i)
	domain = np.arange(-10, 10, 0.01)

	plt.plot(domain, [g(x) for x in domain])
	plt.show()