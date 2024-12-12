""" Main code for the Korteveg-de Vries equation with split-step Fourier method. """

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy as sp
from tqdm import tqdm

# Parametres
L = 50  # Length of the domain
N = 256  # Number of points
dx = L/N
x = np.linspace(0, L, N)

dt = 0.0004  # Time step
t_max = 140
steps = int(t_max / dt)

# Initial conditions
c1 = 0.75 
c2 = 0.4
a1 = 0.33
a2 = 0.65

def u_0(x):
    """ Returns the initial function. """
    return (c1/2)*np.cosh((np.sqrt(c1)/2)*(x-a1*L))**(-2) + (c2/2)*np.cosh((np.sqrt(c2)/2)*(x-a2*L))**(-2)

u0 = u_0(x)

def solution(u_0):
    u_history = np.zeros((steps,N))
    u_history[0] = u_0
    k = np.fft.fftfreq(N)*N
    k3 = 1j*(k*2*np.pi/L)**3
    print()
    with tqdm(total=steps) as pbar:
        for i in range(steps):
            # Partie linéaire
            u_hat = np.fft.fft(u_0)
            u_hat = np.exp(k3 * dt) * u_hat
            u = np.fft.ifft(u_hat).real
            
            # Partie non-linéaire
            u_sq_der = np.fft.ifft(1j*k*np.fft.fft(u**2)*2*np.pi/L)
            
            u = u - 3*dt*u_sq_der
            u = u.real
            u_0 = u
            
            u_history[i] = u
            pbar.update(1)
    print()
    return u_history
 
u_history = solution(u0)

# Analytical solution to compare the numerical solution to
def analytical_sol(t,c,a):
    x = np.linspace(0,L,N)
    u = np.zeros(N)
    for p in range(N):
        u[p] = (np.cosh(np.sqrt(c)*(x[p]-c*t-a*L)/2)**(-2))*c/2
    return u

### Simple plot for 4 cases ###
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

### ----------------------------------------------------------------
mask_x = np.linspace(0, N, N, dtype=int, endpoint=False)
t_plot = np.linspace(0.0, t_max, steps, endpoint=False)
[xx, tt] = np.meshgrid(x[mask_x], t_plot)

fig = plt.figure()
gs = gridspec.GridSpec(3, 4, width_ratios=[1.45, 0.1, 0.20, 1.0])

# Contour plot
ax0 = plt.subplot(gs[:, 0])
contour = ax0.contourf(
    xx, tt, u_history, np.linspace(-0.005, 0.4, 100), cmap='Spectral_r')
ax0.set_title("$N_x = %d, \Delta x = %.4f$" %(256, L/N))
ax0.set_xlabel("$x$")
ax0.set_ylabel("$t$")

# Colorbar
cax = plt.subplot(gs[:, 1])
cbar = fig.colorbar(contour, cax=cax)
cbar.set_ticks(np.linspace(0, 0.3, 6))
plt.show()