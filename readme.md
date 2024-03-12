# Rotational Models for Classical T Tauri Stars
[![DOI](https://zenodo.org/badge/714534647.svg)](https://zenodo.org/doi/10.5281/zenodo.10807253)

This Python code simulates the rotational evolution of young stellar  objects, considering the influence of accretion, magnetic braking via connection star-disk, and accretion-powered stellar winds ([Matt et al., 2012](https://iopscience.iop.org/article/10.1088/0004-637X/745/1/101), [Pinzón et al., 2021](https://iopscience.iop.org/article/10.3847/1538-3881/ac04ae)). The model is based on the evolution of the inertia moment from [Baraffe et al., 2015](https://doi.org/10.1051/0004-6361/201425481) and offers numerical solutions for the system of differential equations governing the star's rotation rate over time.

**Key Features:**

- Provides outputs including time evolution of rotational velocity, and period. Also, it is possible to get torque components by making a modification in the current code.

**Usage:**

1. Specify the stellar parameters, including mass ($`M_{\ast}`$), initial rotation period ($`P^{in}_{rot}`$), initial accretion rate ($`\dot{M}^{in}_{acc}`$), disk lifetime, magnetic field strength ($`B_{\ast}`$), and branching ratio ($`\chi`$ or APSW).
4. Run the `Rotational_models_CTTS` function to obtain time-dependent rotational evolution data.
5. Visualize the rotational evolutionary track using the arrays returned by the function. For enhanced analysis, consider incorporating measurements of $`v\sin(i)`$ and/or $`P_{rot}`$ for individuals or groups of stars alongside the age.

**Dependencies:**

- `numpy`
- `matplotlib`
- `scipy`

**Example Usage:**

Start by extracting the contents of this repository. Then, open an IPython console within the decompressed directory and execute the below code. Alternatively, you can create a new .py file in the repository folder, copy the code below, and run it.

~~~python
import Rotational_models_CTTS
import matplotlib.pyplot as plt

# Set stellar parameters

Mass = 0.5  # Stellar mass in solar masses (Allowed values 0.3, 0.4, 0.5, and so on until, 1.2)
Prot_in = 8.0  # Initial rotation period in days
Macc_in = 1e-8  # Initial accretion rate in solar masses per year (Suggested values between 1e-10 and 1e-6)
Tdisk = 1e7  # Disk lifetime in years (Time when simulation stops)
Bfield = 500  # Magnetic field strength in Gauss (Suggested values between 100 and 3500 G)
betta = 0.01 # Magnetic field coupling parameter (Similar to Matt et al., 2012)
gamma = 1.0 # Magnetic field opening parameter (Similar to Matt et al., 2012)
APSW = 0.01 # Branching ratio parameter (Suggested values between 0.01 and 0.6)

# Run the rotational model

time, vsini, period = Rotational_models_CTTS.Run(Mass, Prot_in, Macc_in, Tdisk, Bfield, betta, gamma, APSW)

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

The model provides an evolutive track of the stellar rotation ($`v\sin(i)`$ and/or $`P_{rot}`$) for a specific input value [$`M_{\ast}`$, $`B_{\ast}`$, $`\chi`$, $`P^{in}_{rot}`$, $`\dot{M}^{in}_{acc}`$]. These output values are particularly valuable as they can be directly compared with observational data (see some applications below).

We provide some outcomes from the use of these models:

In our article (Serna et al., 2023, Submitted to ApJ), we show evolutive tracks for specific cases varying one input parameter while others remain fixed. 

<img src="https://github.com/javiserna/Rotational-models-of-CTTS/blob/main/Plots/Figure_7_A.png?raw=true" width="400"/> <img src="https://github.com/javiserna/Rotational-models-of-CTTS/blob/main/Plots/Figure_7_B.png?raw=true" width="400"/>
<img src="https://github.com/javiserna/Rotational-models-of-CTTS/blob/main/Plots/Figure_7_C.png?raw=true" width="400"/> <img src="https://github.com/javiserna/Rotational-models-of-CTTS/blob/main/Plots/Figure_7_D.png?raw=true" width="400"/>
<img src="https://github.com/javiserna/Rotational-models-of-CTTS/blob/main/Plots/Figure_7_E.png?raw=true" width="400"/>

We also build a grid of models using a parameters space of expected values in CTTS (~$`2.2 \times 10^{6}`$ evolutive tracks in total) and present a comparison of the model's grid with observations.

The grid of models is illustrated as follows as a density plot of $`v\sin(i)=\frac{\pi}{4}v_{rot}`$ as a function of the age, in different values of mass. The black line shows the grid's median $`v\sin(i)`$, while the gray line contains 90\% of the models below it. The color scale represents the number of models per hex bin pixel. The white dots represent measurements of $`v\sin(i)`$ and age for a sample of CTTS from [Hernández et al., 2014](https://iopscience.iop.org/article/10.1088/0004-637X/794/1/36) and [Briceño et al., 2019](https://iopscience.iop.org/article/10.3847/1538-3881/aaf79b).

<img src="https://github.com/javiserna/Rotational-models-of-CTTS/blob/main/Plots/Figure_8_A.png?raw=true" width="400"/> <img src="https://github.com/javiserna/Rotational-models-of-CTTS/blob/main/Plots/Figure_8_B.png?raw=true" width="400"/>
<img src="https://github.com/javiserna/Rotational-models-of-CTTS/blob/main/Plots/Figure_8_C.png?raw=true" width="400"/> <img src="https://github.com/javiserna/Rotational-models-of-CTTS/blob/main/Plots/Figure_8_D.png?raw=true" width="400"/>

The plots show a decreasing trend of the median of $`v\sin(i)`$ with age. Despite the wide range of initial conditions considered for CTTS, most solutions lie on the slow rotator regime with rotation below 10\% of the break-up limit as expected for CTTS. In general, our models agree with observations. For the first time, these models allow us to study how the magnetic field strength (B∗) and branching ratio (χ) change about rotation rates and ages.

**Future Work:**
- As new period measurements in CTTS are available, will be possible to replicate the last plots using evolutionary tracks of $`P_{rot}`$.
  
Measurements of periods in CTTS can be challenging to obtain because of complex brightness variations in the light curves. Despite the challenges, in some cases is possible to get a signature of the rotation in the light curve of these objects.

**References:**

Serna. J, Pinzón. G, Hernández. J et al., 2024 (Accepted for Publication in ApJ). Rotational Evolution of Classical T Tauri Stars: Models and Observations.
