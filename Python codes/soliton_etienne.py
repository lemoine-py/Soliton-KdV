
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Paramètres
L = 50  # Taille du domaine
N = 256  # Nombre de points
dx = L/N
x = np.linspace(0, L, N)
dt = 0.0004  # Pas de temps
t_max = 140  # Temps maximum
steps = int(t_max / dt)

# Conditions initiales
c1 = 0.75 
c2 = 0.4
a1 = 0.33
a2 = 0.65

def u_0(x):
    return (c1/2)*np.cosh((np.sqrt(c1)/2)*(x-a1*L))**(-2) + (c2/2)*np.cosh((np.sqrt(c2)/2)*(x-a2*L))**(-2)

u0 = u_0(x)


def solution(u_0):
    u_history = np.zeros((steps,N))
    u_history[0] = u_0
    k = np.fft.fftfreq(N)*N
    k3 = 1j*(k*2*np.pi/L)**3
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
    return u_history
    
u_history = solution(u0)
plt.plot(x,u_history[-1])

def analytical_sol(t, c, a):
    x = np.linspace(0, L, N)
    u = np.zeros(N)
    for p in range(N):
        # Ajuste x[p] pour être périodique dans [0, L]
        xp = (x[p] - c * t) % L
        u[p] = (np.cosh(np.sqrt(c) * (xp - a * L) / 2) ** -2) * c / 2
    return u


mask_x = np.linspace(0, N, N, dtype=int, endpoint=False)
t_plot = np.linspace(0.0, t_max, steps, endpoint=False)
[xx, tt] = np.meshgrid(x[mask_x], t_plot)


fig = plt.figure(figsize=(15, 7))
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

from matplotlib.animation import FuncAnimation

fig_anim, ax_anim = plt.subplots()

line, = ax_anim.plot([], [])

# Init
ax_anim.set_xlim(x[0], x[-1])
ax_anim.set_ylim(-0.01, 0.4)

ax_anim.set_xlabel("x")
ax_anim.set_ylabel("u")
# ax_anim.set_title(f"Wave Animation, $\lambda={lmb}$")


def animate(i):
    line.set_data(x, u_history[500*i, :])
    return line,

# Create the animation
anim = FuncAnimation(fig_anim, animate, frames=steps//500, blit=True, repeat=True);

# Save the animation as an MP4 file
anim.save('soliton.gif', writer='pillow', fps=25)
