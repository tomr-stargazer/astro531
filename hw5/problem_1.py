"""
Code for Astro 531 HW #5, which involves stars.

"""

from __future__ import division

import numpy as np

import astropy.constants as c
import astropy.units as u



def problem_1a():

	X = 0.707
	Y = 0.274
	mu = 2 / (3*X + Y/2 + 1)

	rho = 1.4 * u.g / u.cm**3
	c_P = 5/2 * c.k_B / (mu*c.u)
	r = c.R_sun / 2
	M = c.M_sun / 2
	T = 10**7 * u.K

	delta_nabla_T = (c.L_sun *  (T/(M*c.G))**(1/2) /
		(4*np.pi * (r**3)/100 * rho * c_P) )**(2/3)

	return delta_nabla_T.decompose() 

