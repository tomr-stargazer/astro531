"""
Problem 1 from Nuria's homework.

"""

from __future__ import division

import numpy as np
import astropy.constants as c
import astropy.units as u

M = c.M_sun
R_0 = 4 * c.R_sun
T_ef = 4000 * u.K

X = 0.707 # hydrogen mass abundance
Y = 0.274 # helium
mu = 2 / (3*X + Y/2 + 1) # mean molecular weight assuming full ionization

def hayashi_radius(time, initial_radius=R_0, effective_temperature=T_ef, mass=M):
	""" Computes the contracted radius after a certain time interval. """

	t = u.Quantity(time, u.s)
	R_0 = initial_radius
	T_ef = effective_temperature
	M = mass

	R_h = (12*np.pi * c.sigma_sb * T_ef**4 * t / (c.G * M**2) + R_0**(-3))**(-1/3)

	return R_h.decompose()

def central_temperature(radius, mean_molecular_weight=mu, mass=M, n=3/2):
	""" Computes the central temperature for a polytrope of index n. """

	R = radius
	M = mass

	# from Clayton Table 2-5, p. 160: n=1.5, column 3 / column 2
	lane_emden_xi_dphi_dxi = 2.71406 / 3.65375

	T_c = mu * u.u / (c.k_B * (n+1) * lane_emden_xi_dphi_dxi) * c.G * M / R

	return T_c.decompose()



