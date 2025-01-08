""" 
___ SOLITON and KORTEWEG-DE VRIES EQUATION ___

Main code to solve numerically the Korteweg-de Vries equation with split-step Fourier method.
This code produces two plots and one animation file.

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy as sp # not used
from tqdm import tqdm # For the progress bar

# Parametres
L = 50  # Length of the domain
N = 256  # Number of points
dx = L/N
x = np.linspace(0, L, N)

dt = 0.0004  # Time step
t_max = 200
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
    """ Returns the solution of the KdV equation. """
    u_history = np.zeros((steps,N))
    u_history[0] = u_0
    k = np.fft.fftfreq(N)*N
    k3 = 1j*(k*2*np.pi/L)**3
    print()
    with tqdm(total=steps) as pbar: # Progress bar
        for i in range(steps):
            # Linear part
            u_hat = np.fft.fft(u_0)
            u_hat = np.exp(k3 * dt) * u_hat
            u = np.fft.ifft(u_hat).real
            
            # Non-linear part
            u_sq_der = np.fft.ifft(1j*k*np.fft.fft(u**2)*2*np.pi/L)
            
            # Euler forward
            u = u - 3*dt*u_sq_der
            u = u.real
            u_0 = u
            
            u_history[i] = u
            pbar.update(1) # Updates the progress bar
    print()
    return u_history

u_history = solution(u0) # Storing the solution array

def analytical_sol(t, c, a):
    x = np.linspace(0, L, N)
    u = np.zeros((N,t_max)) # normally steps not t_max but too heavy
    u_p = np.zeros(N)
    for i in range(t_max): # normally steps not t_max but too heavy
        for p in range(N):
            # Adjust x[i] to be periodic in [0, L]
            xp = (x[p] - c * t) % L
            u_p[p] = (np.cosh(np.sqrt(c) * (xp - a * L) / 2) ** -2) * c / 2
        u[:,i] = u_p[p]
    return u

u_a = analytical_sol(t_max, c1, a1) + analytical_sol(t_max, c2, a2)

### Simple plot for 4 frames

fig, ax = plt.subplots(1, 4, figsize=(20, 5))

ax[0].plot(x, u_history[0], color = "purple")
ax[0].set_title('Initial condition')
ax[0].set_xlabel('$x$')
ax[0].set_ylabel('$u$')
ax[0].set_ylim(-0.01, 0.6)
ax[0].grid()

ax[1].plot(x, u_history[2*steps//20], color = "blue")
ax[1].set_title('t = {:.2f}'.format(2*t_max / 20))
ax[1].set_xlabel('$x$')
ax[1].set_ylim(-0.01, 0.6)
ax[1].grid()

ax[2].plot(x, u_history[3*steps//20], color = "green")
ax[2].set_title('t = {:.2f}'.format(3*t_max / 20))
ax[2].set_xlabel('$x$')
ax[2].set_ylim(-0.01, 0.6)
ax[2].grid()

ax[3].plot(x, u_history[4*steps//20], color = "limegreen")
ax[3].set_title('t = {:.2f}'.format(4*t_max / 20))
ax[3].set_xlabel('$x$')
ax[3].set_ylim(-0.01, 0.6)
ax[3].grid()

plt.tight_layout()
plt.savefig('first_frames_soliton.png')

### Creating the color map
mask_x = np.linspace(0, N, N, dtype=int, endpoint=False)
t_plot = np.linspace(0.0, t_max, steps, endpoint=False)
[xx, tt] = np.meshgrid(x[mask_x], t_plot)

fig = plt.figure(figsize=(15, 7))
gs = gridspec.GridSpec(3, 4, width_ratios=[1.45, 0.1, 0.20, 1.0])

ax0 = plt.subplot(gs[:, 0])
contour = ax0.contourf(
    xx, tt, u_history, np.linspace(-0.005, 0.4, 100), cmap='jet')
ax0.set_title("Korteweg-de Vries equation - two solitons")
y_ticks = [0, 25, 50, 75, 100, 125, 150, 175, 200]
ax0.set_yticks(y_ticks)
ax0.set_xlabel("$x$")
ax0.set_ylabel("$t$")

# Colorbar
cax = plt.subplot(gs[:, 1])
cbar = fig.colorbar(contour, cax=cax)
cbar.set_ticks(np.linspace(0, np.max(u_history), 10))

plt.savefig('cmap_soliton.png')

plt.show()

### Animation
from matplotlib.animation import FuncAnimation

fig_anim, ax_anim = plt.subplots()

line, = ax_anim.plot([], [])

# Init
ax_anim.set_xlim(x[0], x[-1])
ax_anim.set_ylim(-0.01, 0.4)

ax_anim.set_xlabel("x")
ax_anim.set_ylabel("u")
ax_anim.set_title(f"KdV solitons")


def animate(i):
    line.set_data(x, u_history[500*i, :])
    return line,

# Create the animation
anim = FuncAnimation(fig_anim, animate, frames=steps//500, blit=True, repeat=True);

# Save the animation as an MP4 file
anim.save('soliton.gif', writer='pillow', fps=25)