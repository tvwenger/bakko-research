"""
CurrentHIICode.py
Ryan Bakko - 5/13/2024

This program calculates 
    the radio continuum temperature (T_b?),  
    the radio recombination brightness temperature (T_L),
    the optical depth (tau),
    the free-free opacity (also tau?),
    and line center opacity (tau_L)
""" 
import astropy as np
# downloaded: did it work?

# User input
#n_e = float(input("Enter electron density (cm^-3): "))
T_e = float(input("Enter electron temperature (K): "))
nu = float(input("Enter frequency (GHz): "))
EM = float(input("Enter emission measure (pc cm^-6): "))
#weird_nu = float(input("Enter weird frequency (kHz): "))

#tau
def continuum_optical_depth(EM, nu, T_e):
    return((3.014*10^-2)((T_e)^(-3/2))(nu^-2))
    """
    The continuum optical depth (tau), is found using the emission measure, frequency (GHz), and electron temp (T_e = T).
    The optical depth (tau) is simlar to the free-free opactiy (tau). 
    Optical depth requires the free-free gaunt factor, free-free opactiy does not.

    Inputs:
        EM :: Emission measure (pc cm-6)
        nu :: Observation frequency (GHz)
        T_e :: electron temperature (K)

    Outputs:
        tau :: continuum optical depth
    This can be used to calculate T_b
    """
#tau_L
def line_center_opacity(EM, weird_nu, T_e):
    return((1.92*10^3)(T_e^(-5/2))*EM*weird_nu^-1)
    """
    the line center opacity (tau_L) is calculated using emission measure, electronm temp, and a werid frequency value I dont understand.

    Inputs:
        T_e :: electron temp
        EM :: Emission measure (pc cm-6)
        weird_nu :: weird frequency (delta_v/KHz)^-1?
    Outputs:
        tau_L :: line center opacity
    This can be used to calculate T_L
    """
# Function to calculate continuum temperature T_b, the brightness temperature of free-free emission
def continuum_temperature(T_e,tau):
    return (T_e(1-e^-tau))
    """
    T_b, one of the main goals of this program

    Inputs:
        T_e :: electron temp
        tau :: optical depth
    Outputs:
        T_b :: continuum temperature, the brightness temperature of free-free emission
    """

# Function to calculate radio recombination brightness temperature T_L
def radio_recombination_brightness(T_e,tau_L):
   return(T_e(1-e^-tau_L))
   """
    T_L is the other main calculation of this program.
    T_L = T_e * tau_L 
    OR T_L = T_e(1-e^-(tau_L)) ?

    Inputs:
        T_e :: electron temp
        tau_L :: line center opacity
    Outputs:
        T_L :: brightness temperature of recombination emission line at its center frequency 
    """

# Calculating and printing results
print("Continuum temperature: {:.2f} K".format(continuum_temperature(T_e, continuum_optical_depth)))
print("Radio recombination brightness: {:.2f} K".format(radio_recombination_brightness(T_e, line_center_opacity)))
