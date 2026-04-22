import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import pickle

dim = 4
k = 0.1

def gaussian(sigma):
    return lambda x : -1 / (np.sqrt(2 * np.pi) * sigma) * np.exp(-0.5 * (x / sigma) ** 2)

V = gaussian(0.2)

def func(x, y):
    der = ((-1j) ** dim) * (k ** dim - V(x)) * y[0]
    sol = np.append(y[1:dim], der)
    return sol

roots = np.array([1j * k * np.exp(2 * np.pi * 1j * l / dim) for l in range(dim)], dtype=np.cdouble)
print("roots", roots)

x0 = -10
xf = 10
coefficients = []

for momentum in roots:
	iv = np.array([momentum ** i * np.exp(momentum * x0) for i in range(dim)], dtype=np.cdouble)
     
	sol = sc.integrate.solve_ivp(func, [x0, xf], iv, t_eval=np.linspace(x0, xf, 300))

	a = np.full((dim, dim), 1j, dtype=np.cdouble)
	for der in range(dim):
		for (i, mom) in enumerate(roots):
			a[der][i] = mom ** der * np.exp(mom * xf)
	b = sol.y.T[-1]
	print("a and b")
	print(a.shape)
	print(b.shape)
	sol = np.linalg.solve(a, b)

	coefficients.append(sol)

	if False:
		wavefunction = [sol.y[0, i] for i in range(sol.y.shape[1])]
		print(len(wavefunction))
		plt.scatter(sol.t, np.real(wavefunction))
		plt.plot(sol.t, np.real(np.exp(momentum * sol.t)), color='orange')
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