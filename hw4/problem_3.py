"""

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import trapz

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

# equation 5
eta_r = solar_model['L/Lsun'] / solar_model['M/Msun']

P_r = solar_model['P']

chi_R_eta_bar_r = np.zeros_like(chi_R_r)

for i in range(len(r)):

    # integrate from the current radius OUTWARDS
    chi_eta_integrand = (chi_R_r * eta_r)[i:]

    chi_R_eta_bar_r[i] = 1 / P_r[i] * -trapz(chi_eta_integrand, P_r[i:])

print chi_R_eta_bar_r[0]
print chi_R_eta_bar_r[-1]

# make a beta

thing = c.L_sun / (4 * np.pi * c.c * c.G * c.M_sun) * chi_R_eta_bar_r

print thing.decompose().unit

beta_r = 1 - c.L_sun / (4 * np.pi * c.c * c.G * c.M_sun) * chi_R_eta_bar_r

# problem4


anisotropic_term = (3 / (4 * np.pi) * solar_model['L/Lsun']*c.L_sun /
                    (4*np.pi * (solar_model['R/Rsun'] * c.R_sun)**2) )

anisotropic_term_surface = (3 / (4 * np.pi) * c.L_sun /
                    (4*np.pi * (solar_model['R/Rsun'] * c.R_sun)**2) )

isotropic_term = (c.c / (4 * np.pi) * (1 - beta_r) * solar_model['P']*u.dyn/u.cm**2)

def problem4():
    """ 
    Make some plots.

    """

    fig = plt.figure()

    plt.plot(r, (anisotropic_term/isotropic_term).decompose(), 
             label='L = local luminosity')
    plt.plot(r, (anisotropic_term_surface/isotropic_term).decompose(), 
             label='L = surface luminosity')

    plt.semilogy()
    plt.title("Problem 4: Ratio of anisotropic to isotropic radiation in the sun")

    plt.xlabel(r"r/R$_\odot$")
    plt.ylabel("anisotropic / isotropic")

    plt.legend(loc='lower right')
             
    return fig
