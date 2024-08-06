import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as c

# Forward model function
def forward_model(physical_params):
    T_e, nu_0, EM, non_thermal_fwhm = physical_params

    def fwhm(T_e, nu_0):
        return ((8 * np.log(2) * c.k_B) / c.c**2)**0.5 * (T_e / c.m_p)**0.5 * nu_0

    def continuum_optical_depth(EM, nu, T_e):
        return 3.28e-7 * (T_e.to("K").value / 1.0e4) ** -1.35 * (nu.to("GHz").value ** -2.1) * EM.to("cm-6 pc").value

    def line_center_opacity(EM, fwhm, T_e):
        return 1.92e3 * (T_e.to("K").value ** -2.5) * (EM.to("cm-6 pc").value) / (fwhm.to("GHz").value)
    
    def line_profile(nu, nu_0, tau_l, fwhm_freq):
         return tau_l * np.exp(-4.0 * np.log(2.0) * (nu - nu_0)**2.0 / fwhm_freq**2.0)

    def brightness_temperature(tau_c, T_e):
        return T_e * (1 - np.exp(-tau_c))

    T_e = T_e * u.K
    nu_0 = nu_0 * u.GHz
    EM = EM * u.pc / u.cm**6
    non_thermal_fwhm = non_thermal_fwhm * u.km / u.s

    fwhm_val = fwhm(T_e, nu_0)
    tau_c = continuum_optical_depth(EM, nu_0, T_e)
    tau_l = line_center_opacity(EM, fwhm_val, T_e)
    brightness_temp = brightness_temperature(tau_c, T_e)

    return brightness_temp.value

# Adding noise function
def add_noise(spectrum, loc=0.0, scale=1.0):
    noise = np.random.normal(loc=loc, scale=scale, size=spectrum.shape)
    return spectrum + noise
    #noise_spec = np.random.normal(loc=0.0, scale=noise, size=len(nu)) * u.K

# Loss function
def loss_function(physical_params, observed_spectrum):
    predicted_spectrum = forward_model(physical_params)
    residuals = observed_spectrum - predicted_spectrum
    return residuals

# Prepare your data (generate synthetic observations)
freqs = np.linspace(1.0, 10.0, 100) * u.GHz  # Frequency range for the spectrum - New
true_physical_params = [8000.0, 7.0, 1000.0, 25.0]  # Replace with true physical parameters
true_spectrum = forward_model(true_physical_params)
observed_spectrum = add_noise(true_spectrum, loc=0.0, scale=0.1)  # Adjust scale


# Set up initial parameters
initial_guess = [7500.0, 6.5, 950.0, 20.0]  # Initial guess for the physical parameters

# Run least squares fitting
result = least_squares(loss_function, initial_guess, args=(observed_spectrum,))

# Extract and use the results
estimated_physical_params = result.x
print("Estimated Physical Parameters:", estimated_physical_params)

# Plot the observed and fitted spectrum
plt.plot(true_spectrum, label='True Spectrum', linestyle='--')
plt.plot(observed_spectrum, label='Observed Spectrum', linestyle=':')
plt.plot(forward_model(estimated_physical_params), label='Fitted Spectrum', linestyle='-')
plt.legend()
plt.xlabel('Frequency (GHz)')
plt.ylabel('Brightness Temperature (K)')
plt.title('Observed vs Fitted Spectrum')
plt.show()
