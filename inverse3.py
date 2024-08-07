"""
inverse3.py
Experiments with forward modeling (despite the name!)
Ryan Bakko & Trey Wenger - August 2024

TODO: add more complete docstrings to every function!
"""

import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as c
import CurrentHIICode

def forward_model(freqs, physical_params, nu_0):
    """
    Predict an observed radio recombination line and radio continuum
    spectrum.

    Inputs:
        freqs :: 1-D array of scalars (GHz)
            Frequencies at which to evaluate the spectrum
        physical_params :: 4-length array of scalars
            physical_params[0] = Electron temperature (K)
            physical_params[1] = Emission measure (pc cm-6)
            physical_params[2] = Non-thermal FWHM line width in velocity units (km s-1)
        nu_0 :: scalar (GHz)
            RRL rest frequency

    Returns:
        spectrum :: 1-D array of scalars (K)
            Brightness temperature spectrum
    """
    T_e, EM, non_thermal_fwhm = physical_params

    # add units
    T_e = T_e * u.K
    nu_0 = nu_0 * u.GHz
    EM = EM * u.pc / u.cm**6
    non_thermal_fwhm = non_thermal_fwhm * u.km / u.s

    # evaluate physics
    non_thermal_fwhm_freq = nu_0 * non_thermal_fwhm / c.c
    thermal_fwhm_freq = CurrentHIICode.fwhm(T_e, nu_0)
    fwhm_freq = np.sqrt(thermal_fwhm_freq**2.0 + non_thermal_fwhm_freq**2.0)
    tau_c = CurrentHIICode.continuum_optical_depth(EM, freqs, T_e)
    tau_l = CurrentHIICode.line_center_opacity(EM, fwhm_freq, T_e)
    tau_l_profile = CurrentHIICode.line_profile(freqs, nu_0, tau_l, fwhm_freq)
    brightness_temp = CurrentHIICode.brightness_temperature(T_e, tau_c + tau_l_profile)

    return brightness_temp.to('K').value

# Adding noise function
def add_noise(spectrum, loc=0.0, scale=1.0):
    noise = np.random.normal(loc=loc, scale=scale, size=spectrum.shape)
    return spectrum + noise

# Loss function
def loss_function(physical_params, freqs, nu_0, observed_spectrum):
    predicted_spectrum = forward_model(freqs, physical_params, nu_0)
    residuals = observed_spectrum - predicted_spectrum
    return residuals

# Prepare your data (generate synthetic observations)
freqs = np.linspace(6.99, 7.01, 1000) * u.GHz  # Frequency range for the spectrum - New
true_physical_params = [8000.0, 1000.0, 25.0]  # Replace with true physical parameters
nu_0 = 7.0 # RRL rest frequency (GHz)
true_spectrum = forward_model(freqs, true_physical_params, nu_0)
noise = 0.001
observed_spectrum = add_noise(true_spectrum, loc=0.0, scale=noise)  # Adjust scale

# Set up initial parameters
initial_guess = [7500.0, 950.0, 20.0]  # Initial guess for the physical parameters

# Run least squares fitting
result = least_squares(loss_function, initial_guess, args=(freqs, nu_0, observed_spectrum,))

# Extract and use the results
estimated_physical_params = result.x
print("Estimated Physical Parameters:", estimated_physical_params)

# Plot the observed and fitted spectrum
plt.plot(observed_spectrum, 'k-', label='Observed Spectrum')
plt.plot(true_spectrum, 'r-', label='True Spectrum')
plt.plot(forward_model(freqs, estimated_physical_params, nu_0), 'b-', label='Fitted Spectrum')
plt.legend()
plt.xlabel('Frequency (GHz)')
plt.ylabel('Brightness Temperature (K)')
plt.title('Observed vs Fitted Spectrum')
plt.show()
