"""
CurrentHIICode2.py

Simulate radio continuum and radio recombination line observations of
Galactic HII regions.

Ryan Bakko - 5/20/2024
Ryan Bakko & Trey Wenger - July 2024
"""

import sys

from astropy import units as u
from astropy import constants as c
import numpy as np

def fwhm(T_e):
    """
    Find the spectral width using Temperature
    Equation 7.35 from Essential Radio Astronomy textbook

    Inputs:
        T_e :: electron temp (K)
        M :: mass (g)
        k :: Boltzmann constant? (erg K**-1)
        c :: spped of light (cm/s)
        v_0 :: inital (GHz)

    Output:
        delta_v :: fwhm (Hz)
    """ 
    #struggling with formatting this equation into my code
    return (
        ((8 * ln(2) * k) / c**2 )**0.5
        * (T_e.to("K").value / M)**0.5
        * v_0
    )

def continuum_optical_depth(EM, nu, T_e):
    """
    Calculate the radio continuum optical depth.
    Equation 4.60 from Condon & Ransom ERA textbook.

    Inputs:
        EM :: Emission measure (with units)
        nu :: Observation frequency (with units)
        T_e :: electron temperature (with units)

    Outputs:
        tau :: continuum optical depth
    """
    return (
        3.28e-7
        * ((T_e.to("K").value / 1.0e4) ** -1.35)
        * (nu.to("GHz").value ** -2.1)
        * EM.to("cm-6 pc").value
    )


def line_center_opacity(EM, fwhm, T_e):
    """
    Calculate the radio recombination line optical depth at the line center.
    Equation 7.96 from Condon & Ransom ERA textbook.

    Inputs:
        EM :: Emission measure (with units)
        fwhm :: FWHM line width in frequency units (with units)
        T_e :: electron temperature (with units)

    Outputs:
        tau_L :: radio recombination line optical depth at the line center
    """
    return (
        (1.92e3)
        * (T_e.to("K").value ** -2.5)
        * (EM.to("cm-6 pc").value)
        * (fwhm.to("kHz").value ** -1)
    )


def brightness_temperature(T_e, tau):
    """
    Calculates the brightness temperature.

    Inputs:
        Te :: electron temperature (with units)
        tau :: Optical depth (unitless)

    Outputs:
        T_b :: brightness temperature (with units)
    """
    return T_e * (1 - np.exp(-tau))


def main(T_e, nu, EM, fwhm):
    """
    Simulate radio continuum and radio recombination line observations of a Galactic HII region.

    Inputs:
        T_e :: electron temperature (with units)
        nu :: frequency (with units)
        EM :: emission measure (with units)
        fwhm :: FWHM line width in velocity units (with units)

    Outputs:
        Nothing
    """
    # TODO: calculate line width from temperature- above
    print(f"T_e = {T_e.to('K'):.1f}")
    print(f"nu = {nu.to('MHz'):.1f}")
    print(f"EM = {EM.to('pc cm-6'):.1e}")
    print(f"fwhm = {fwhm.to('km/s'):.1f}")
    print()

# Iterate over the array of frequencies- NEW
    for freq in nu:
        print(f"nu = {freq.to('MHz'):.1f}")

    # convert FWHM to frequency units
    fwhm_freq = nu * fwhm / c.c

    # Continuum optical depth
    tau_c = continuum_optical_depth(EM, nu, T_e)
    print(f"tau_c = {tau_c:.3e}")

    # RRL line center optical depth
    tau_l = line_center_opacity(EM, fwhm_freq, T_e)
    print(f"tau_l = {tau_l:.3e}")

    # Brightness temperature
    TB_c = brightness_temperature(T_e, tau_c)
    TB_l = brightness_temperature(T_e, tau_l)
    print(f"TB_c = {TB_c.to('mK'):.1f}")
    print(f"TB_l = {TB_l.to('mK'):.1f}")

    # line-to-continuum ratio
    line_to_cont = TB_l / TB_c
    print(f"line-to-continuum ratio = {line_to_cont:.3f}")


if __name__ == "__main__":
    # Read inputs from command line
    # e.g. python CurrentHIICode.py 8000.0 7.0 1000.0 25.0
    if len(sys.argv) > 1:
        T_e = float(sys.argv[1]) * u.K
        nu = float(sys.argv[2]) * u.GHz
        EM = float(sys.argv[3]) * u.pc / u.cm**6
        fwhm = float(sys.argv[4]) * u.km / u.s
    else:
        print("Using default parameters")
        T_e = 8000.0 * u.K
        nu = 7.0 * u.GHz
        EM = 1000.0 * u.pc / u.cm**6
        fwhm = 25.0 * u.km / u.s

    main(T_e, nu, EM, fwhm)