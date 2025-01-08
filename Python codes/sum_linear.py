"""
Sum of the two solutions of each term in the Korteweg-de Vries equation.
Demonstrates that the solution for the KdV equation is not linear.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy as sp # not used
from tqdm import tqdm # For the progress bar
from matplotlib.animation import FuncAnimation

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

def analytical_sol(t, c, a, steps):
    x = np.linspace(0, L, N)
    u = np.zeros((N, steps)) # normally steps not t_max but too heavy
    
    for i in range(steps): # normally steps not t_max but too heavy
        u_p = np.zeros(N)
        for p in range(N):
            # Adjust x[i] to be periodic in [0, L]
            xp = (x[p] - c * i) % L
            u_p[p] = (np.cosh(np.sqrt(c) * (xp - a * L) / 2) ** -2) * c / 2
        u[:,i] = u_p
    return u

def u_array(u, steps):
    u_array = []
    for i in range(steps):
        u_array.append(u[:,i])
    return u_array

### Simple plot for 4 frames

def four_frames(u,u_max,t_end,steps):
    fig, ax = plt.subplots(1, 4, figsize=(20, 5))

    ax[0].plot(x, u[:,0], color = "purple")
    ax[0].set_title('Initial condition')
    ax[0].set_xlabel('$x$')
    ax[0].set_ylabel('$u$')
    ax[0].set_ylim(-0.01, u_max)
    ax[0].grid()

    ax[1].plot(x, u[:,2*t_end//20], color = "blue")
    ax[1].set_title('t = {:.2f}'.format(2*t_end / 20))
    ax[1].set_xlabel('$x$')
    ax[1].set_ylim(-0.01, u_max)
    ax[1].grid()

    ax[2].plot(x, u[:,3*t_end//20], color = "green")
    ax[2].set_title('t = {:.2f}'.format(3*t_end / 20))
    ax[2].set_xlabel('$x$')
    ax[2].set_ylim(-0.01, u_max)
    ax[2].grid()

    ax[3].plot(x, u[:,4*t_end//20], color = "limegreen")
    ax[3].set_title('t = {:.2f}'.format(4*t_end / 20))
    ax[3].set_xlabel('$x$')
    ax[3].set_ylim(-0.01, u_max)
    ax[3].grid()

    plt.tight_layout()
    plt.savefig('first_frames_linear.png')
    plt.show()

### Creating the color map

def colormap(t_end, steps, u_array, u_max):
    mask_x = np.linspace(0, N, N, dtype=int, endpoint=False)
    t_plot = np.linspace(0.0, t_end, steps, endpoint=False)
    [xx, tt] = np.meshgrid(x[mask_x], t_plot)

    fig = plt.figure(figsize=(15, 7))
    gs = gridspec.GridSpec(3, 4, width_ratios=[1.45, 0.1, 0.20, 1.0])


    ax0 = plt.subplot(gs[:, 0])
    contour = ax0.contourf(
        xx, tt, u_array, np.linspace(-0.005, u_max, 100), cmap='jet')
    ax0.set_title("KdV equation - Sum of linear solutions")
    y_ticks = [0, 25, 50, 75, 100, 125, 150, 175, 200]
    ax0.set_yticks(y_ticks)
    ax0.set_xlabel("$x$")
    ax0.set_ylabel("$t$")

    # Colorbar
    cax = plt.subplot(gs[:, 1])
    cbar = fig.colorbar(contour, cax=cax)
    cbar.set_ticks(np.linspace(0, np.max(u_array), 10))

    plt.savefig('cmap_soliton.png')
    plt.show()

### Animation

def gif_creator(u_array, x, steps, u_max):
    """Create a gif file with the animation of the solution of the KdV equation"""
    fig_anim, ax_anim = plt.subplots()

    line, = ax_anim.plot([], [])

    # Init
    ax_anim.set_xlim(x[0], x[-1])
    ax_anim.set_ylim(-0.01, u_max)

    ax_anim.set_xlabel("x")
    ax_anim.set_ylabel("u")
    ax_anim.set_title(f"KdV solitons - Sum of linear solutions")

    def animate():
        line.set_data(x, u_array)
        return line,

    # Create the animation
    anim = FuncAnimation(fig_anim, animate, frames=steps, blit=True, repeat=True);

    # Save the animation as a gif file
    anim.save('soliton.gif', writer='pillow', fps=25)


### Calling the functions

u_a = analytical_sol( t_max, c1, a1, t_max) + analytical_sol(t_max, c2, a2, t_max)

u_sol = u_array(u_a, t_max)

four_frames(u_a, 0.6, t_max, steps)

colormap(t_max, t_max, u_sol, 0.8)

gif_creator(u_sol, x, steps, 0.8)