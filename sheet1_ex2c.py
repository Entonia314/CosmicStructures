import numpy as np
import matplotlib.pyplot as plt

H0 = 100
c = 299792458

# Function to calculate scale factor a(eta) in open universe
def a_open(eta, Omega_m):
    return Omega_m / (2 * (1 - Omega_m)) * (np.cosh(np.sqrt(1 - Omega_m) * eta) - 1)


# Function to calculate time t(eta) in open universe
def t_open(eta, Omega_m):
    numerator = np.sinh(np.sqrt(1-Omega_m) * eta)
    denominator = np.sqrt(1-Omega_m)
    return Omega_m / (2 * (1-Omega_m)) * (numerator / denominator - eta)


# Function to calculate scale factor a(eta) in closed universe
def a_closed(eta, Omega_m):
    return Omega_m / (2 * (Omega_m - 1)) * (1 - np.cos(eta))


# Function to calculate time t(eta) in closed universe
def t_closed(eta, Omega_m):
    return Omega_m * (eta - np.sin(eta)) / ((Omega_m - 1)**(3/2)*2*H0)


# Function to calculate scale factor a(eta) in flat universe
def a_flat(eta, Omega_m):
    A = Omega_m * H0**2 / (2 * c**2)
    return A / 2 * eta**2


# Function to calculate time t(eta) in flat universe
def t_flat(eta, Omega_m):
    A = Omega_m * H0**2 / (2 * c**2)
    return A / 6*c * eta**3


# Generate values of conformal time eta
eta_values = np.linspace(0, 10, 100)

# Define values of Omega_m to plot
Omega_m_values = [0.3, 1.00001, 1.3]

# Plot scale factor a(eta) for each Omega_m
plt.figure(figsize=(10, 6))
"""for Omega_m in Omega_m_values:
    print("Omega_m: ", Omega_m)
    print(Omega_m/(2*(1-Omega_m)))
    print(np.sqrt(1-Omega_m))
    a_values = a_open(eta_values, Omega_m)
    plt.plot(eta_values, a_values, label=f'$\\Omega_m = {Omega_m}$')"""
plt.plot(eta_values, a_open(eta_values, 0.3), label=f'$\\Omega_m = {0.3}$')
plt.plot(eta_values, a_flat(eta_values, 1), label=f'$\\Omega_m = {1}$')
plt.plot(eta_values, a_closed(eta_values, 1.3), label=f'$\\Omega_m = {1.3}$')

plt.title('Scale Factor $a(\\eta)$ vs. Conformal Time $\\eta$')
plt.xlabel('Conformal Time $\\eta$')
plt.ylabel('Scale Factor $a(\\eta)$')
plt.legend()
plt.grid(True)
plt.show()
