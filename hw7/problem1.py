"""
Problem 1 from Nuria's homework.

"""

from future import __division__

import numpy as np
import astropy.constants as c
import astropy.units as u

M = c.M_sun
R_0 = 4 * c.R_sun
T_ef = 4000 * u.K


def hayashi_radius(initial_radius, time, effective_temperature, mass):
	""" Computes the contracted radius after a certain time interval. """

	R_0 = initial_radius
	t = u.Quantity(time, u.s)
	T_ef = effective_temperature
	M = mass

	R_h = (12 * np.pi * c.sigma_sb * T_ef**4 * t / (c.G * M**2) + R_0**(-1/3))**(-1/3)

	return R_h

def central_temperature(mass, radius):
	""" Computes the central temperature for a polytrope of index n. """

	pass
