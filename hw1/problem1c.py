"""
Code to make calculations & plots for Problem 1c.

"Plot the spectral energy distribution 
log (nu F_nu) vs. log (lambda).
Also for comparison overplot 
log (nu F_nu) vs. log (lambda)
calculated with fluxes not corrected for reddening."

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from astropy import constants
import astropy.units as u

wavelength_array = np.array([0.36, 0.44, 
                             0.55, 0.64, 
                             0.79, 1.25, 
                             1.65, 2.2])*u.um

frequency_array = (constants.c / wavelength_array).to('Hz')

dereddened_magnitude_array = np.array([13.85, 12.86,
                                       11.62, 10.82,
                                       10.03, 9.01,
                                       8.31, 8.20])

observed_magnitude_array = np.array([14.61, 13.50,
                                     12.11, 11.23,
                                     10.33, 9.15,
                                     8.40, 8.25])

reddening_array = observed_magnitude_array - dereddened_magnitude_array

dereddened_flux_array = np.array(
    [5.22e-29, 3.06e-28, 8.19e-28, 1.44e-27,
     2.48e-27, 4.02e-27, 4.97e-27, 3.44e-27] ) * u.W / (u.m)**2 / u.Hz

observed_flux_array = dereddened_flux_array / 10**(reddening_array / 2.5)

def plot_SED():

    plt.plot(wavelength_array, frequency_array*dereddened_flux_array,
             'bo--')

    plt.plot(wavelength_array, frequency_array*observed_flux_array, 
             'ro:')

    plt.loglog()
    plt.xlim(0.35, 2.3)
    plt.ylim(1e-14, 2e-12)

    plt.xlabel(r"Wavelength ($\mu m$)")
    plt.ylabel(r"$\nu F_{\nu}$ (W m$^{-2}$)")

    plt.xticks(wavelength_array, 
               ('U', 'B', 'V', 'R', 'I', 'J', 'H', 'K'))

    plt.text(0.4, 6.5e-13, "Dereddened flux", color='b')
    plt.text(1, 6.5e-13, "Observed flux", color='r')

    plt.title("SED of V826 Tau")
