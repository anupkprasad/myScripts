import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import k, m_u

def maxwell_boltzmann_speed(v, T, m):
    """Computes the Maxwell-Boltzmann speed distribution."""
    factor = (m / (2 * np.pi * k * T)) ** (3 / 2)
    return 4 * np.pi * v**2 * factor * np.exp(-m * v**2 / (2 * k * T))

# Constants
mass_N2 = 28.0 * m_u  # Nitrogen molecule mass in kg
T_values = [100, 300, 1000]  # Different temperatures in Kelvin

# Speed range
v = np.linspace(0, 2000, 500)

# Plot
plt.figure(figsize=(8, 6))
for T in T_values:
    plt.plot(v, maxwell_boltzmann_speed(v, T, mass_N2), label=f"T = {T} K")

plt.xlabel("Speed (m/s)")
plt.ylabel("Probability Density")
plt.title("Maxwell-Boltzmann Speed Distribution")
plt.legend()
plt.grid()
plt.show()

