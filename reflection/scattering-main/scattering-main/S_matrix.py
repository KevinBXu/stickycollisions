import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import lu
from scipy import special
import sys

orig_stdout = sys.stdout
f = open('log.log', 'w')
sys.stdout = f

np.set_printoptions(suppress=True)
np.set_printoptions(linewidth=np.inf)
np.set_printoptions(precision=3)
import pickle


with open('S-matrix.pkl', 'rb') as file:
    S = pickle.load(file)

even = S[2, 0] + S[2, 1]

decay1 = even / -S[2, 2]

odd = S[2, 0] - S[2, 1]

decay2 = odd / -S[2, 2]

even_wavefunction = np.array([1, 1, decay1, 0])
odd_wavefunction = np.array([1, -1, decay2, 0])

x1 = S @ even_wavefunction 
x2 = S @ odd_wavefunction
print(x1)
print(x2)

left = np.array([[x1[0], x2[0]], [1, -1]])
right = np.array([[1, 1], [x1[1], x2[1]]])
S2 = left @ np.linalg.inv(right)
print("S2")
print(S2)
print(S2 @ [1, 1])
eigenvalues, eigenvectors = np.linalg.eig(S2)
print("Eigenvalues")
print(eigenvalues)
print("Eigenvectors")
print(eigenvectors)
inverse = np.array([[1, 1], [1, -1]])


# a = S[0, 0]
# b = S[0, 1]
# c = S[1, 0]
# d = S[1, 1]

# c1 = a - b * c / d
# c2 = b / d
# c3 = -c / d
# c4 = 1 / d

# S3 = np.array([[c1, c2], [c3, c4]])
# print("S3")
# print(S3)
# print(S3 @ [1, 1])
# eigenvalues, eigenvectors = np.linalg.eig(S3)
# print("Eigenvalues")
# print(eigenvalues)
# print("Eigenvectors")
# print(eigenvectors)


f.close()
exit()

c1 = 1
c2 = 1
c3 = -(S @ np.array([c1, c2, 0, 0]))[2] / S[2, 2]
input = np.array([1, 1, c3, 0])
output = S @ input
print(output)
x1 = output[:2]
y1 = [c1, c2]

c1 = 1
c2 = -1
c3 = -(S @ np.array([c1, c2, 0, 0]))[2] / S[2, 2]
input = np.array([1, -1, c3, 0])
output = S @ input
print(output)
x2 = output[:2]
y2 = [c1, c2]

left = np.array([[x1[0], x2[0]], [1, -1]])
right = np.array([[1, 1], [x1[1], x2[1]]])
S = left @ np.linalg.inv(right)
inverse = np.array([[1, 1], [1, -1]])

eigenvalues, eigenvectors = np.linalg.eig(S)