""" Main code for the Korteveg-de Vries equation with split-step Fourier method. """

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from tqdm import tqdm

# Parametres
L = 50  # Length of the domain
N = 256  # Number of points
dx = L / N
x = np.linspace(0, L, N, endpoint=False)

dt = 0.0004  # Time step
t_max = 140
steps = int(t_max / dt)

# Initial conditions
c1, c2 = 0.75, 0.4
a1, a2 = 0.33, 0.65
u = (c1/2)*np.cosh((np.sqrt(c1)/2)*(x-a1*L))**(-2) + (c2/2)*np.cosh((np.sqrt(c2)/2)*(x-a2*L))**(-2)

# Fourier vectors
k = 2 * np.pi * np.fft.fftfreq(N, d=dx)
k3 = 1j * k**3

print()
# Resolution with split-step Fourier
u_history = [u.copy()]
with tqdm(total=steps) as pbar:
    for _ in range(steps):
        # linear part in spectral space
        u_hat = np.fft.fft(u)
        u_hat = np.exp(k3 * dt) * u_hat
        u = np.fft.ifft(u_hat).real
        
        # non-linear part in physical space
        u_sq_der = np.fft.ifft(1j*k*np.fft.fft(u**2))
    
        u = u - 3*dt*u_sq_der
        u = u.real
    
        # storage for visualization
        u_history.append(u.copy())
        pbar.update(1)
print()

fig, ax = plt.subplots(1, 4, figsize=(20, 5))
ax[0].plot(x, u_history[0])
ax[0].set_title('Initial condition')
ax[0].set_xlabel('x')
ax[0].set_ylabel('u')
ax[0].grid()

ax[1].plot(x, u_history[steps//4])
ax[1].set_title('t = {:.2f}'.format(t_max/4))
ax[1].set_xlabel('x')
ax[1].set_ylabel('u')
ax[1].grid()

ax[2].plot(x, u_history[steps//2])
ax[2].set_title('t = {:.2f}'.format(t_max/2))
ax[2].set_xlabel('x')
ax[2].set_ylabel('u')
ax[2].grid()

ax[3].plot(x, u_history[-1])
ax[3].set_title('t = {:.2f}'.format(t_max))
ax[3].set_xlabel('x')
ax[3].set_ylabel('u')
ax[3].grid()
plt.show()
