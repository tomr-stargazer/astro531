"""

"""

from __future__ import division

import numpy as np

import astropy
import astropy.units as u
import astropy.constants as c

from astro531.hw4.problem_2 import solar_model

# let's calculate (chi_R eta) bar

r = solar_model['R/Rsun']
X_r = solar_model['X']
Y_r = solar_model['Y(He4)'] + solar_model['He3']
rho_r = solar_model['Rho'] * u.g / u.cm**3
T_r = solar_model['T'] * u.K

# see equation 6 in the pset
chi_R_r = (4e4 * u.cm**2/u.g * (1 - X_r - Y_r) * (1 + X_r) * 
           (rho_r / (u.g / u.cm**3)) / (T_r/(1e6*u.K))**3.5).to('cm2 g-1')

print chi_R_r[0]

def chi_R_eta_bar(r):
    pass
    
