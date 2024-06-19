
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

from astropy.constants import c
from astropy import units as u
import numpy as np



#tau
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
    # TODO: organize function comments like above -- Done
    # 3.014e-2
    # products: not (a)(b) is: a * b
    # (EM / cm-6 pc) TODO: fix units everywhere -- Done
    # Temperature xample: T.to("K").value/1.0e4
    return((3.014e-2)*((T_e)**(-3/2))*(nu**-2)*EM.to("EM / cm-6 pc").value)


#tau_L, radio recombination line optical depth at the line center
def line_center_opacity(EM, weird_nu, T_e):
    
    """
    The radio recombination line optical depth at the line center (tau_L) is calculated using emission measure, 
    electronm temp, and werid nu (the width of the spectral line).

    Inputs:
        T_e :: electron temp (K)
        EM :: Emission measure (EM / cm-6 pc)
        weird_nu :: width of the spectral line (delta_v/KHz)**-1

    Outputs:
        tau_L :: radio recombination line optical depth at the line center
    This can be used to calculate T_L
    """
    return((1.92e3)*(T_e**(-5/2))*EM.to("EM / cm-6 pc")*weird_nu**-1)


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
    # User input
    #n_e = float(input("Enter electron density (cm**-3): "))
    T_e = float(input("Enter electron temperature (K): ")) * u.K # TODO: do this elsewhere ---DONE!
    nu = float(input("Enter frequency (GHz): ")) * u.GHz
    EM = float(input("Enter emission measure (EM / cm-6 pc): ")) * (u.EM /(u.cm**-6 * u.pc))
    weird_nu = float(input("Enter width of the spectral line (delta_v / KHz)**-1: ")) * (u.delta_v / u.KHz**-1)

    # Convert input values to astropy units
    T_e = T_e * u.K
    nu = nu * u.GHz
    EM = EM * (u.pc * u.cm**-6)
    weird_nu = weird_nu * u.kHz

    # Calculating and printing results
    tau_c = continuum_optical_depth(EM, nu, T_e)
    tau_l = line_center_opacity(EM, weird_nu, T_e)
    print("Continuum temperature: {:.2f} K".format(brightness_temperature(T_e, tau_c)))
    print("Radio recombination brightness: {:.2f} K".format(brightness_temperature(T_e, tau_l)))
    #print(f"Continuum temperature: {T_b:.2f}")
    #print(f"Radio recombination brightness: {T_L:.2f}")

if __name__ == "__main__":
    main()
