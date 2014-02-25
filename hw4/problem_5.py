"""
Problem 5 of Nuria's HW4. Adiabatic exponent gamma3 stuff.

"""

from __future__ import division

import numpy as np

import astropy
import astropy.units as u
import astropy.constants as c

chi_H = 13.6 * u.eV

def ionization_fraction_y(temperature, density):
    """
    Ionization fraction 'y' following Clayton's notation.

    """

    T = u.quantity.Quantity(temperature, u.K)
    rho = u.quantity.Quantity(density, u.g / u.cm**3)

    saha_g = ((2 * np.pi * c.m_e * c.k_B * T)**(3/2) / c.h**3 * np.exp(- chi_H / (c.k_B * T) )).decompose().value

    Hplus = 1/2 * ( ( ((c.u + c.m_e)*saha_g)**2 + 4*saha_g*rho)**(1/2) -
                    (c.u + c.m_e)*saha_g )

    H_plus_Hplus = rho / c.u

    y = Hplus / H_plus_Hplus

    return y


def gamma3_minus_one(temperature, density ):
    """

    """

    T = u.quantity.Quantity(temperature, u.K)

    y = ionization_fraction_y(temperature, density)

    D = y*(1-y)/( (2-y) * (1+y) )

    threehalfs_plus_chiH_over_kT = 3/2 + chi_H/(c.k_B * T) 

    reduced_gamma = ((2 + 2*D * threehalfs_plus_chiH_over_kT ) /
                     (3 + 2*D * t`xhreehalfs_plus_chiH_over_kT**2 ) ) 

    return reduced_gamma

