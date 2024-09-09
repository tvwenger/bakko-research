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

import matplotlib.pyplot as plt

def fwhm(T_e, nu_0):
    """
    Calculate the FWHM of the spectral line
    Find the spectral width using Temperature
    Equation 7.35 from Essential Radio Astronomy textbook

    Inputs:
        T_e :: electron temp (w units)
        nu_0 :: rest freq (w units)

    Output:
        delta_nu :: fwhm (w units)
    """ 
   
    return (
        ((8 * np.log(2) * c.k_B) / c.c**2 )**0.5
        * (T_e / c.m_p)**0.5
        * nu_0
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

def line_profile(nu, nu_0, tau_l, fwhm_freq):
    """
    evaluate the gaussian line profile
    
    Inputs:
        nu :: frequencies (with units)
        nu_0 :: rest frequency (with units)
        tau_l :: line-center opacity
        fwhm_freq :: FWHM in frequency units (with units)

    Outputs:
        line_prof :: line profile evaluated at nu
    """
    return tau_l * np.exp(-4.0 * np.log(2.0) * (nu - nu_0)**2.0 / fwhm_freq**2.0)

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


def main(T_e, nu_0, EM, non_thermal_fwhm):
    """
    Simulate radio continuum and radio recombination line observations of a Galactic HII region.

    Inputs:
        T_e :: electron temperature (with units)
        nu_0 :: rest frequency (with units)
        EM :: emission measure (with units)
        non_thermal_fwhm :: Non thermal FWHM line width in velocity units (with units)

    Outputs:
        Nothing
    """
    # Input parameters
    print(f"T_e = {T_e.to('K'):.1f}")
    print(f"nu_0 = {nu_0.to('MHz'):.1f}")
    print(f"EM = {EM.to('pc cm-6'):.1e}")
    print(f"non_thermal_fwhm = {non_thermal_fwhm.to('km/s'):.1f}")
    print()

    # Calculate the thermal FWHM and Iterate over the array of frequencies
    thermal_fwhm = fwhm(T_e, nu_0)

    #Convert Non-Thermal Line width from km/s to Hz
    non_thermal_fwhm_freq = nu_0 * non_thermal_fwhm / c.c

    #Total line width calculation using both thermal_fwhm and non_thermal_fwhm
    fwhm_freq = np.sqrt(thermal_fwhm**2+non_thermal_fwhm_freq**2)

    #Range of frequencies for the profile
    nu = np.linspace(nu_0 - 5 * fwhm_freq, nu_0 + 5 * fwhm_freq, 1000 )

    # Continuum optical depth
    tau_c_spec = continuum_optical_depth(EM, nu, T_e)
    #print(f"tau_c = {tau_c:.3e}")

    # RRL line center optical depth
    tau_l = line_center_opacity(EM, fwhm_freq, T_e)
    # print(f"tau_l = {tau_l:.3e}")

    # Optical depth spectrum for the line profile
    tau_l_spec = line_profile(nu, nu_0, tau_l, fwhm_freq)

    # Brightness temperature spectrum for continuum and line
    TB_c_spec = brightness_temperature(T_e, tau_c_spec)
    TB_l_spec = brightness_temperature(T_e, tau_l_spec)
    # print(f"TB_c = {TB_c.to('mK'):.1f}")
    # print(f"TB_l = {TB_l.to('mK'):.1f}")

    # line-to-continuum ratio
    # line_to_cont = TB_l / TB_c
    # print(f"line-to-continuum ratio = {line_to_cont:.3f}")

    # Plot the results - Having trouble with line labels
    fig, ax = plt.subplots()
    ax.plot(nu, TB_c_spec + TB_l_spec, 'k-', label='Total Brightness Temp')
    ax.plot(nu, TB_c_spec, 'r--', label='Continuum brightness Temp')

    # Labels and units for axes
    ax.set_xlabel(f'Frequency (GHz)')
    ax.set_ylabel(f'Brightness Temperature (K)')
    
    # Legend
    # ax.legend(loc='upper right', fontsize=12, frameon=True, shadow=True)

    # Title
    ax.set_title('Radio Continuum and Recombination Line Observations')

    # Display the plot
    plt.show(block=True)

    

if __name__ == "__main__":
    # Read inputs from command line
    # e.g. python CurrentHIICode.py 8000.0 7.0 1000.0 25.0
    if len(sys.argv) > 1:
        T_e = float(sys.argv[1]) * u.K
        nu_0 = float(sys.argv[2]) * u.GHz
        EM = float(sys.argv[3]) * u.pc / u.cm**6
        non_thermal_fwhm = float(sys.argv[4]) * u.km / u.s

    else:
        print("Using default parameters")
        T_e = 8000.0 * u.K
        nu_0 = 7.0 * u.GHz
        EM = 1000.0 * u.pc / u.cm**6
        non_thermal_fwhm = 25.0 * u.km / u.s

    main(T_e, nu_0, EM, non_thermal_fwhm)