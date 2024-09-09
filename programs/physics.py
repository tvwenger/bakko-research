"""
physics.py
Physical utility functions.
Trey Wenger - September 2024
"""

import pytensor.tensor as pt
from astropy import constants as c

_C = c.c.to("km/s")
_K_B = c.k_B.to("km2 s-2 g K-1")
_M_P = c.m_p.to("g")


def calc_thermal_fwhm(electron_temp, rest_freq):
    """
    Calculate the thermal FWHM of the spectral line
    Equation 7.35 from Essential Radio Astronomy textbook

    Inputs:
        electron_temp :: electron temp (K)
        rest_freq :: rest freq (MHz)

    Output:
        delta_nu :: fwhm (MHz)
    """
    return pt.sqrt(8.0 * pt.log(2.0) * _K_B / _C**2.0) * pt.sqrt(electron_temp / _M_P) * rest_freq


def calc_continuum_optical_depth(freqs, electron_temp, emission_measure):
    """
    Calculate the radio continuum optical depth.
    Equation 4.60 from Condon & Ransom ERA textbook.

    Inputs:
        freqs :: Observed frequencies (MHz)
        electron_temp :: electron temperature (K)
        emission_measure :: emission measure (pc cm-6)

    Outputs:
        continuum_optical_depth :: continuum optical depth
    """
    return 3.28e-7 * ((electron_temp / 1e4) ** -1.35) * (freqs / 1000.0) ** -2.1 * emission_measure


def calc_line_center_optical_depth(emission_measure, fwhm_freq, electron_temp):
    """
    Calculate the radio recombination line optical depth at the line center.
    Equation 7.96 from Condon & Ransom ERA textbook.

    Inputs:
        emission_measure :: Emission measure (pc cm-6)
        fwhm_freq :: FWHM line width in frequency units (MHz)
        electron_temp :: electron temperature (K)

    Outputs:
        line_center_optical_depth :: radio recombination line optical depth at the line center
    """
    return 1.92e3 * (electron_temp**-2.5) * emission_measure * (1000.0 * fwhm_freq) ** -1.0


def calc_optical_depth_line_profile(freqs, rest_freq, line_center_optical_depth, fwhm_freq):
    """
    Evaluate the optical depth line profile

    Inputs:
        freqs :: frequencies (MHz)
        rest_freq :: rest frequency (MHz)
        line_center_optical_depth :: line-center opacity
        fwhm_freq :: FWHM in frequency units (MHz)

    Outputs:
        line_prof :: optical line profile evaluated at freqs
    """
    return line_center_optical_depth * pt.exp(-4.0 * pt.log(2.0) * (freqs - rest_freq) ** 2.0 / fwhm_freq**2.0)


def calc_brightness_temperature(electron_temp, optical_depth):
    """
    Calculates the brightness temperature.

    Inputs:
        electron_temp :: electron temperature (K)
        optical_depth :: Optical depth (unitless)

    Outputs:
        T_b :: brightness temperature (K)
    """
    return electron_temp * (1.0 - pt.exp(-optical_depth))
