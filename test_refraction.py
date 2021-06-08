# -*- coding: utf-8 -*-
"""
This is the test file ofr the python script refraction_shift.py. 

Author       : Hanae Labriji
Institute    : Office National d'Etudes et de Recherches Aerospatiales ONERA
Last modified: 2021-06-08
"""

from RefractionShift.refraction_shift import refraction
import unittest
import numpy as np

class TestFunctionRefraction(unittest.TestCase):
    def test_AngularShift(self):
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
        
        ## Check if Ra has the expected value
        self.assertAlmostEqual(Ra, 0.000288, places = 5)

    
    def test_LateralShift(self):
        T_0 = 273.15 # temperature at the observer (K)
        P_0 = 100000 # pressure at the observer (Pa)
        h0 = 0 # altitude at the observer (m)
        lmbda =  550e-9 # Wavelength (m)
        z0_deg = 45 # zenith angle (deg) 
        
        ## Create an instance of the class LateralShift
        Refraction = refraction(T_0, P_0, h0)
        
        ## Compute the refraction angle at 45 (deg)
        z0 = np.deg2rad(z0_deg) # zenith angle (rad)
        LatShift = Refraction.get_LateralShift(z0, lmbda)
        
        ## Check if Ra has the expected value
        self.assertAlmostEqual(LatShift, 3.261232, places = 5)
        
if __name__ == '__main__':
    unittest.main()