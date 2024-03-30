import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.misc import derivative

Omega_r = 0.00008
Omega_m = 0.3
Omega_Lambda = 0.7
H0 = 70

lambda_matter_equ = 0.7646
matter_radiation_equ = 0.00027193


def Hubble_a(a):
    return H0 * np.sqrt(Omega_r * a**(-4) + Omega_m * a**(-3) + Omega_Lambda)

def t_of_a(a):
    res = np.zeros_like(a)
    for i, ai in enumerate(a):
        t, err = quad(lambda ap: 1.0/(ap*Hubble_a(ap)), 0, ai)
        res[i] = t
    return res


def radiation_epoch(a):
    return a**2 / (2 * np.sqrt(Omega_r))


def matter_epoch(a):
    return a**(3/2) / (H0 * np.sqrt(Omega_m)*(3/2)) * H0


def lambda_epoch(a):
    return np.log(a/lambda_matter_equ) / (H0 * np.sqrt(Omega_Lambda)) * H0


a = np.logspace(-8, 1, 100)


plt.plot(radiation_epoch(a), a, label="Radiation Epoch")
plt.plot(matter_epoch(a), a, label="Matter Epoch")
plt.plot(lambda_epoch(a), a, label="Lambda Epoch")
plt.loglog(t_of_a(a) * H0, a, label="Numerical")
plt.vlines([lambda_matter_equ, matter_radiation_equ], 0, 10, colors=["pink", "purple"])
plt.xlabel('t * H0')
plt.ylabel('a(t)')
plt.title('Numerical Solution vs analytical asymptotic solutions')
plt.legend()
plt.show()

