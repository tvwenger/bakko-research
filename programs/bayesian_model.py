"""
bayesian_model.py
Demonstration of a Bayesian model
Trey Wenger - September 2024
"""

import pymc as pm
import pytensor.tensor as pt
import arviz as az
import numpy as np

import physics


def predict_spectrum(freqs, rest_freq, electron_temp, emission_measure, nonthermal_fwhm):
    """
    Predict the RRL + continuum brightness temperature spectrum.

    Inputs:
        freqs :: 1-D array of scalars (MHz)
            Observed frequencies
        rest_freq :: scalar (MHz)
            Rest frequency
        electron_temp :: scalar (K)
            Electron temperature
        emission_measure :: scalar (pc cm-6)
            Emission measure
        nonthermal_fwhm :: scalar (km/s)
            Non-thermal FWHM line width

    Returns:
        spectrum :: 1-D array of scalars (mK)
            Brightness temperature spectrum
    """
    # Non-thermal FWHM in frequency units (MHz)
    nonthermal_fwhm_freq = rest_freq * nonthermal_fwhm / physics._C

    # Thermal FWHM in frequency units (MHz)
    thermal_fwhm_freq = physics.calc_thermal_fwhm(electron_temp, rest_freq)

    # Total line width (MHz)
    fwhm_freq = pt.sqrt(thermal_fwhm_freq**2.0 + nonthermal_fwhm_freq**2.0)

    # Continuum optical depth
    tau_c = physics.calc_continuum_optical_depth(freqs, electron_temp, emission_measure)

    # Line center RRL optical depth
    tau_l_center = physics.calc_line_center_optical_depth(emission_measure, fwhm_freq, electron_temp)

    # RRL optical depth line profile
    tau_l = physics.calc_optical_depth_line_profile(freqs, rest_freq, tau_l_center, fwhm_freq)

    # Brightness temperature (K)
    brightness_temp = physics.calc_brightness_temperature(electron_temp, tau_c + tau_l)
    return 1000.0 * brightness_temp


def main():
    # Generate synthetic data
    freqs = np.linspace(6990.0, 7010.0, 1000)  # MHz
    rest_freq = 7000.0  # MHz
    true_electron_temp = 8000.0  # K
    true_emission_measure = 1000.0  # pc cm-6
    true_nonthermal_fwhm = 25.0  # km s-1
    true_noise = 1.0  # mK
    # Predicted brightness temperature spectrum (mK)
    true_spectrum = predict_spectrum(
        freqs, rest_freq, true_electron_temp, true_emission_measure, true_nonthermal_fwhm
    ).eval()
    # Add noise (mK)
    obs_spectrum = true_spectrum + true_noise * np.random.randn(*true_spectrum.shape)

    # Define model
    with pm.Model() as model:
        # Prior distribution on model parameters
        log10_electron_temp = pm.Normal("log10_electron_temp", mu=4.0, sigma=1.0)  # K
        log10_emission_measure = pm.Normal("log10_emission_measure", mu=5.0, sigma=1.0)  # pc cm-6
        log10_nonthermal_fwhm = pm.Normal("log10_nonthermal_fwhm", mu=1.0, sigma=1.0)  # km s-1
        noise = pm.HalfNormal("noise", sigma=10.0)  # mK

        # Evaluate the model
        predicted_brightness_temp = predict_spectrum(
            freqs, rest_freq, 10.0**log10_electron_temp, 10.0**log10_emission_measure, 10.0**log10_nonthermal_fwhm
        )

        # Likelihood function
        _ = pm.Normal("brightness_temp", mu=predicted_brightness_temp, sigma=noise, observed=obs_spectrum)

    # Sample model with MCMC
    with model:
        trace = pm.sample()

    print(az.summary(trace))


if __name__ == "__main__":
    main()
