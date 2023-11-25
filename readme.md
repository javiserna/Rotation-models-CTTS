# Rotational Models for Classical T Tauri Stars

This Python code simulates the rotational evolution of young stellar  objects, considering the influence of accretion, magnetic braking via the presence of a protoplanetary disk, and accretion-powered stellar winds ([Matt et al., 2012](https://iopscience.iop.org/article/10.1088/0004-637X/745/1/101), [Pinzón et al., 2021](https://iopscience.iop.org/article/10.3847/1538-3881/ac04ae)). The model is based on the evolution of the inertia moment from [Baraffe et al. 2015](https://doi.org/10.1051/0004-6361/201425481) and offers numerical solutions for the system of differential equations governing the star's rotation rate over time.

**Key Features:**

- Provides outputs including time evolution of rotational velocity, period, torque components, and stellar radii.

**Usage:**

1. Specify the stellar parameters, including mass, rotation period, accretion rate, disk lifetime, and magnetic field strength.
2. Run the `Rotational_model` function to obtain time-dependent rotational evolution data.
3. Analyze and visualize the results using the returned arrays.

**Dependencies:**

- `numpy`
- `matplotlib`
- `scipy`

**Example Usage:**

Begin by decompressing this repository. Next, open an IPython console and execute the provided code directly or create a new .py file, copy the content, and run it within the main repository directory.

~~~python
import Rotational_models_CTTS
import matplotlib.pyplot as plt

# Set stellar parameters

Mass = 1.0  # Stellar mass in solar masses (Allowed values 0.3, 0.4, 0.5, and so on until, 1.2)
Prot = 5.0  # Initial rotation period in days
Macc = 1e-7  # Accretion rate in solar masses per year (Suggested values between 1e-10 and 1e-6)
Tdisk = 1e7  # Disk lifetime in years (Time when simulation stops)
Bfield = 2500  # Magnetic field strength in Gauss (Suggested values between 100 and 3500 G)
betta = 0.01 # Magnetic field coupling parameter (Similar to Matt et al., 2012)
gamma = 1.0 # Magnetic field opening parameter (Similar to Matt et al., 2012)
APSW = 0.1 # Branching ratio parameter (Suggested values between 0 and 0.6)

# Run the rotational model

time, vsini, period = Rotational_models_CTTS.Run(Mass, Prot, Macc, Tdisk, Bfield, betta, gamma, APSW)

# Plot the results
plt.figure(1)
plt.plot(time/1e6, vsini)
plt.xlabel('Time (Myr)')
plt.ylabel('Rotational Velocity (km/s)')

plt.figure(2)
plt.plot(time/1e6, period)
plt.xlabel('Time (Myr)')
plt.ylabel('Rotation Period (days)')
plt.show()
~~~

The simulation outputs provide the evolutive track of the stellar rotation (v sin(i) and/or Prot). These values are particularly valuable as they can be directly compared with observational data.

**References:**

Serna. J, Pinzón. G, Hernández. J et al., 2023 (Submitted to ApJ). Rotational Evolution of Classical T Tauri Stars: Models and Observations.
