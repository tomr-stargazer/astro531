"""
Problem 5 of Nuria's HW4. Adiabatic exponent gamma3 stuff.

"""

from __future__ import division

import numpy as np

import astropy
import astropy.units as u
import astropy.constants as c

def ionization_fraction_y():
    """
    Ionization fraction 'y' following Clayton's notation.

    """

    # y = H+ / (H + H+)

    return y


def gamma3_minus_one(temperature, density ):
    """

    """

    T = temperature
    chi_H = 13.6 * u.eV

    D = y*(1-y)/( (2-y) * (1+y) )

    threehalfs_plus_chiH_over_kT = 3/2 + chi_H/(c.k_B * T) 

    reduced_gamma = (2 + 2*D * threehalfs_plus_chiH_over_kT ) / (3 + 2*D * threehalfs_plus_chiH_over_kT**2 ) 

    return reduced_gamma

