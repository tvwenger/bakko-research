
"""
CurrentHIICode.py
Ryan Bakko - 5/20/2024

This program calculates 
    the radio continuum temperature (T_b),  
    the radio recombination brightness temperature (T_L),
    the optical depth (tau),
    the  free-free continuum optical depth (also tau, explained below),
    and radio recombination line optical depth at the line center (tau_L)
""" 

from astropy import units as u
import numpy as np


#tau, continuum optical depth
def continuum_optical_depth(EM, nu, T_e):
    """
    The continuum optical depth (tau), is found using the emission measure, frequency (GHz), and electron temp (T_e = T).
    The optical depth (tau) is simlar to the  free-free continuum optical depth (tau). 
    Optical depth requires the free-free gaunt factor,  free-free continuum optical depth does not.

    Inputs:
        EM :: Emission measure (pc cm-6)
        nu :: Observation frequency (GHz)
        T_e :: electron temperature (K)

    Outputs:
        tau :: continuum optical depth
    This can be used to calculate T_b
    """
    # Rewrote equation without free free gaunt factor
    return((3.28e-7) * ((T_e.to("K").value)**(-1.35)) * (nu.to("GHz").value**-2.1) * EM.to("cm-6 pc").value)


#tau_L, radio recombination line optical depth at the line center
def line_center_opacity(EM, spectral_line_width, T_e):
    
    """
    The radio recombination line optical depth at the line center (tau_L) is calculated using emission measure, 
    electronm temp, and werid nu (the width of the spectral line).

    Inputs:
        T_e :: electron temp (K)
        EM :: Emission measure (EM / cm-6 pc)
        spectral_line_width :: width of the spectral line (delta_v/KHz)**-1

    Outputs:
        tau_L :: radio recombination line optical depth at the line center
    This can be used to calculate T_L
    """
    return((1.92e3) * ((T_e.to("K").value/1e4)**(-5/2)) * (EM.to("cm-6 pc").value) * ((spectral_line_width.to("kHz").value)**-1))


# Function to calculate either continuum temperature T_b, the brightness temperature of free-free emission OR radio recombination brightness temperature, T_L
def brightness_temperature(T_e, tau):
    """
    Calculates either continuum or recombination line brightness temperature.

    Inputs:
        Te :: electron temperature (K)
        tau :: Optical depth (either continuum or line) 
        OR tau_L :: radio recombination line optical depth at the line center

    Outputs:
        T_b OR T_L :: brightness temperature (K) of either continuum or line
    """
    return T_e * (1 - np.exp(-tau))
    #how to choose? TODO: Give either “tau” or “tau_L” and have one equation to calculate either line or continuum.

def main():
    """
    # User input
    #n_e = float(input("Enter electron density (cm**-3): "))
    T_e = float(input("Enter electron temperature (K): ")) 
    nu = float(input("Enter frequency (GHz): "))
    EM = float(input("Enter emission measure (cm-6 pc): "))
    spectral_line_width = float(input("Enter width of the spectral line (kHz): "))
    """
    T_e = 8000.0
    nu = 9.0
    EM = 1000.0
    spectral_line_width = 1000.0

    # Convert input values to astropy units
    T_e = T_e * u.K
    nu = nu * u.GHz
    EM = EM * (u.pc * u.cm**-6)
    spectral_line_width = spectral_line_width * u.kHz

    print(f"tau_cont = {continuum_optical_depth(EM, nu, T_e):.3e}")

    # Calculating and printing results
    tau_c = continuum_optical_depth(EM, nu, T_e)
    tau_l = line_center_opacity(EM, spectral_line_width, T_e)
    print("Continuum temperature: {:.2f} ".format(brightness_temperature(T_e, tau_c)))
    print("Radio recombination brightness: {:.2f} ".format(brightness_temperature(T_e, tau_l)))
    #print(f"Continuum temperature: {T_b:.2f}")
    #print(f"Radio recombination brightness: {T_L:.2f}")

if __name__ == "__main__":
    main()
