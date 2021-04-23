# -*- coding: utf-8 -*-
"""
This is a minimal working example for the python script refraction_shift.py. 

Author       : Hanae Labriji
Institute    : Office National d'Etudes et de Recherches Aerospatiales ONERA
Last modified: 2021-04-22
"""

from src.refraction_shift import refraction
import numpy as np

##############################################################################
## Define the main constants of the observation
T_0 = 273.15 # temperature at the observer (K)
P_0 = 100000 # pressure at the observer (Pa)
h0 = 0 # altitude at the observer (m)
lmbda =  550e-9 # Wavelength (m)
z0_deg = 45 # zenith angle (deg) 

## Create an instance of the class LateralShift
Refraction = refraction(T_0, P_0, h0)

## Compute the refraction angle at 45 (deg)
z0 = np.deg2rad(z0_deg) # zenith angle (rad)
Ra = Refraction.get_AngularShift(z0, lmbda)
print("The refraction angle at the zenith angle : ", np.rad2deg(z0), 
      "deg is equal to ", round(np.rad2deg(Ra)*3600, 2), "arcseconds.")

## Compute the lateral shift at 45 (deg)
LatShift = Refraction.get_LateralShift(z0, lmbda)
print("The lateral shift is equal to ", round(LatShift, 2), "meters.")

## Compute the first order approximation of the lateral shift
LatShift1 = Refraction.get_approx1_LateralShift(z0,lmbda)
print("First order approximation of the lateral shift is equal to ", 
      round(LatShift1, 2), "meters.")

## Compute the second order approximation of the lateral shift
LatShift2 = Refraction.get_approx2_LateralShift(z0,lmbda)
print("Second order approximation of the lateral shift is equal to ", 
      round(LatShift2, 2), "meters.")

## Check the result for z0 equal to 45 deg
if abs(LatShift - 3.261232) > 1e-6 or abs(Ra - 0.000288) > 1e-6:
    print("Computation Error ! Check your input parameters.")
