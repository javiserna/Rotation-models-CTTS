# Rotational Models for Classical T Tauri Stars

This Python code simulates the rotational evolution of young stellar  objects, considering the influence of accretion, magnetic braking via the presence of a protoplanetary disk, and accretion-powered stellar winds ([Matt et al., 2012](https://iopscience.iop.org/article/10.1088/0004-637X/745/1/101), [Pinzón et al., 2021](https://iopscience.iop.org/article/10.3847/1538-3881/ac04ae)). The model is based on the evolution of the inertia moment from [Baraffe et al. 2015](https://doi.org/10.1051/0004-6361/201425481) and includes numerical solutions for the system of differential equations governing the star's rotation rate over time.

**Key Features:**

- Utilizes Baraffe stellar evolution data for initial conditions.
- Numerically integrates differential equations to model accretion, magnetic braking, and wind effects.
- Adaptive step size adjustment during integration for accuracy and efficiency.
- Provides outputs including time evolution of rotational velocity, torque components, radii, mass, and period.

**Usage:**

1. Specify the stellar parameters, including mass, rotation period, accretion rate, disk lifetime, and magnetic field strength.
2. Run the `Rotational_model` function to obtain time-dependent rotational evolution data.
3. Analyze and visualize the results using the returned arrays.

**Dependencies:**

- `numpy`
- `matplotlib`
- `scipy`

**Example Usage:**

~~~python
from Rotational-Models-CTTS import Rotational_model
import matplotlib.pyplot as plt

# Set stellar parameters

Mass = 1.0  # Stellar mass in solar masses
Prot = 8.0  # Initial rotation period in days
Macc = 1e-7  # Accretion rate in solar masses per year
Tdisk = 1e6  # Disk lifetime in years
Bfield = 2500  # Magnetic field strength in Gauss
betta = 0.01 # Remain fixed
gamma = 1.0 # Remain fixed
APSW = 0.1 # Branching ratio

# Run the rotational model

time, vsini, period = Rotational_model(Mass, Prot, Macc, Tdisk, Bfield, betta, gamma, APSW)

# Plot the results
plt.figure(1)
plt.plot(time, vsini)
plt.xlabel('Time (years)')
plt.ylabel('Rotational Velocity (km/s)')
plt.legend()

plt.figure(2)
plt.plot(time, period)
plt.xlabel('Time (years)')
plt.ylabel('Rotation Period (day)')
plt.show()
~~~

The simulation outputs include the time evolution of the stellar rotational velocity (vsini) and rotational period. These values are particularly valuable as they can be directly compared with observational data.

**References:**

Serna. J, Pinzón. G, Hernández. J et al., 2023 (Submitted to ApJ). Rotational Evolution of Classical T Tauri Stars: Models and Observations.
