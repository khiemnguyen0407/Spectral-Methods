"""
Band-limited interpolation
"""

import numpy as np
import matplotlib.pyplot as plt

h = 1
xmax = 10
x = np.arange(-xmax, xmax+h, h)
xx = np.linspace(-xmax, xmax+h/20, 151)

plt.figure(figsize=(5,9))
for iplt in range(3):
    plt.subplot(3, 1, iplt+1)
    if iplt == 0:
        v = np.array(x == 0, dtype=float)       # delta function
    elif iplt == 1:
        v = np.array(np.abs(x) <= 3, dtype=float)               # square wave
    elif iplt == 2:
        v = np.maximum(0, 1 - np.abs(x)/3)      # hat function
        
    plt.plot(x, v, 'k.', markersize=10)
    p = np.zeros(len(xx))
    for i in range(len(x)):
        p = p + v[i] * np.sin(np.pi * (xx - x[i])/h) / (np.pi * (xx - x[i]) /h)
    plt.plot(xx, p, linewidth=2)
    plt.xticks(x[::2])
    plt.yticks(np.linspace(0, 1, 3))
    plt.grid(True)
    
plt.tight_layout
