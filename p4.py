"""
Periodic spectral differentiation
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import toeplitz

# Set up grid and differentiation matrix
N = 24; h = 2*np.pi/N; x = h*np.arange(1, N+1);
vec = np.arange(1, N)
col = np.append(np.array([0]), 0.5*(-1)**vec / np.tan(0.5*h*vec))
row = np.zeros(col.shape)
row[1:] = col[N-1:0:-1]
D = toeplitz(col, row)

# Differentiation of a hat function
v = np.maximum(0, 1 - abs(x - np.pi)/2)

plt.figure(figsize=(8,6), dpi=80)
plt.subplot(3,2,1)
plt.plot(x, v, '-o', markersize=3); plt.grid(True)
plt.axis([0, 2*np.pi, -0.5, 1.5])
plt.title('function')
plt.subplot(3,2,2)
plt.plot(x, np.dot(D, v), '-o', markersize=3); plt.grid(True)
plt.axis([0, 2*np.pi, -1, 1])
plt.title('spectral derivative')

# Differentiation of exp(sin(x))
v = np.exp(np.sin(x)); vprime = np.cos(x) * v
plt.subplot(3,2,3)
plt.plot(x, v, '-o', markersize=3); plt.grid(True)
plt.subplot(3,2,4)
plt.plot(x, np.dot(D, v), '-o', markersize=3); plt.grid(True)
error = np.linalg.norm(np.dot(D, v) - vprime, np.inf)
plt.axis([0, 2*np.pi, -2, 2])
plt.text(2.2, 1.4, 'max error = ' + str(np.round(error, 18)))