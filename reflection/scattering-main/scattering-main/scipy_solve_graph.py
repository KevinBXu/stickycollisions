import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import pickle

dim = None

def gaussian(sigma):
    return lambda x : -1 / (np.sqrt(2 * np.pi) * sigma) * np.exp(-0.5 * (x / sigma) ** 2)

V = gaussian(0.2)

def func(x, y):
    der = ((-1j) ** dim) * (k ** dim - V(x)) * y[0]
    sol = np.append(y[1:dim], der)
    return sol


x0 = -5
xf = 5
t_range = np.linspace(x0, xf, 300)

def run(k):
	coefficients = []
	roots = np.array([1j * k * np.exp(2 * np.pi * 1j * l / dim) for l in range(dim)], dtype=np.cdouble)
	for momentum in roots:
		iv = np.array([momentum ** i * np.exp(momentum * x0) for i in range(dim)], dtype=np.cdouble)
		
		sol = sc.integrate.solve_ivp(func, [x0, xf], iv, method='BDF')

		a = np.full((dim, dim), 1j, dtype=np.cdouble)
		for der in range(dim):
			for (i, mom) in enumerate(roots):
				a[der][i] = mom ** der * np.exp(mom * xf)
		b = sol.y.T[-1]
		assert(np.abs(sol.t[-1] - xf) < 0.001)
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

	even_output = -compute(np.zeros(dim // 2), True)[1:dim // 2]
	even_sol = np.linalg.solve(coeffs, even_output)
	x1 = compute(even_sol, True)

	odd_output = -compute(np.zeros(dim // 2), False)[1:dim // 2]
	coeffs = S[1:dim // 2, 1:dim // 2]
	odd_sol = np.linalg.solve(coeffs, odd_output)
	x2 = compute(odd_sol, False)

	left = np.array([[x1[0], x2[0]], [1, -1]])
	right = np.array([[1, 1], [x1[dim // 2], x2[dim // 2]]])
	S2 = left @ np.linalg.inv(right)
	# print("S2")
	# print(S2)
	# print(S2 @ [1, 1])
	eigenvalues, eigenvectors = np.linalg.eig(S2)
	eigenvalues = np.conj(eigenvalues)
	print(eigenvalues)
	angle = np.angle(eigenvalues) / np.pi
	for i in range(2):
		if angle[i] < 0:
			angle[i] += 2
	numbers = dim * angle / 2
	# print(angle)
	even_vector = np.array([0.707+0.j, 0.707+0.j])
	# print(eigenvectors[:,0], eigenvectors[:,1])
	# print(numbers)
	if np.linalg.norm(eigenvectors[:,0] - even_vector) < 0.5:
		return numbers[0], numbers[1]
	else:
		return numbers[1], numbers[0]
	

dim = 8
V = gaussian(0.5)
title = "scipyGaussian,d=" + str(dim)
pot = [V(x) for x in t_range]
plt.plot(t_range, pot)
plt.show()

even_eigs = []
odd_eigs = []
ks = []
for n in range(0, 50):
    k = 1 - 0.02 * n
    ks.append(k)
    print(k)
    even_eig, odd_eig = run(k)
    even_eigs.append(even_eig)
    odd_eigs.append(odd_eig)
# for i in range(1, len(even_eigs)):
#     if even_eigs[i] - even_eigs[i - 1] > dim // 2:
#         even_eigs[i] -= dim
#     elif even_eigs[i] - even_eigs[i - 1] < -dim // 2:
#         even_eigs[i] += dim
# for i in range(1, len(odd_eigs)):
#     if odd_eigs[i] - odd_eigs[i - 1] > dim // 2:
#         odd_eigs[i] -= dim
#     elif odd_eigs[i] - odd_eigs[i - 1] < -dim // 2:
#         odd_eigs[i] += dim
plt.plot(ks, even_eigs, label='even')
plt.plot(ks, odd_eigs, label='odd')
plt.xlabel("k")
plt.ylabel("m")
plt.title(title) 
plt.legend()
plt.savefig(f"graphs/{title}.pdf")
plt.show()