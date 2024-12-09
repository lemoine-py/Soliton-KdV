import numpy as np
import matplotlib.pyplot as plt

# Paramètres
L = 50  # Taille du domaine
N = 256  # Nombre de points
dx = L / N
x = np.linspace(0, L, N, endpoint=False)
dt = 0.0004  # Pas de temps
t_max = 140  # Temps maximum
steps = int(t_max / dt)

# Conditions initiales
c1, c2 = 0.75, 0.4
a1, a2 = 0.33, 0.65
u = (c1/2)*np.cosh((np.sqrt(c1)/2)*(x-a1*L))**(-2) + (c2/2)*np.cosh((np.sqrt(c2)/2)*(x-a2*L))**(-2)

# Vecteurs de Fourier
k = 2 * np.pi * np.fft.fftfreq(N, d=dx)
k3 = 1j * k**3

# Résolution avec split-step Fourier
u_history = [u.copy()]
for _ in range(steps):
    # Partie linéaire en espace spectral
    u_hat = np.fft.fft(u)
    u_hat = np.exp(k3 * dt) * u_hat
    u = np.fft.ifft(u_hat).real
    
    # Partie non-linéaire en espace physique
    u_sq_der = np.fft.ifft(1j*k*np.fft.fft(u**2))
    
    u = u - 3*dt*u_sq_der
    u = u.real

    # Stockage pour visualisation
    u_history.append(u.copy())
plt.plot(x,u_history[-1])
