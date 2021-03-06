"""
Code for Astro 531 HW #5, which involves stars.

"""

from __future__ import division

import numpy as np

import astropy.constants as c
import astropy.units as u

# Various constants that are used later in the math.
X = 0.707 # hydrogen mass abundance
Y = 0.274 # helium
mu = 2 / (3*X + Y/2 + 1) # mean molecular weight assuming full ionization

rho = 1.4 * u.g / u.cm**3 # mass density
c_P = 5/2 * c.k_B / (mu*c.u) # heat capacity
r = c.R_sun / 2 # radius under consideration 
M = c.M_sun / 2 # mass inside r
T = 10**7 * u.K # temperature

l = r / 10 # mixing length at r

gamma = 5/3 # adiabatic index
P = rho/(mu * c.u) * c.k_B * T # pressure

ln_Lambda = 10 # cutoff factor in nu (Spitzer 1962)
# Coefficient of kinematical viscosity
nu = (2.21e-15 * u.g / u.cm / u.s) * (T/u.K)**(5/2) / (rho * ln_Lambda)


def problem_1a():
	""" 
	Calculate the excess of the temperature gradient over the adiabatic gradient delta nabla T.

	"""

	delta_nabla_T = (c.L_sun *  (T/(M*c.G))**(1/2) /
		(4*np.pi * (r**3)/100 * rho * c_P) )**(2/3)

	# unit hacking because there was a weird s0 unit term that wouldn't go away
	return u.Quantity(delta_nabla_T.decompose().value, u.K / u.m).to('K/cm')

def	problem_1b():
	"""	Calculate the temperature excess delta nabla T x l. """

	delta_nabla_T = problem_1a()

	excess = delta_nabla_T * l

	return excess.decompose()

def problem_1c():
	""" Calculate the velocity of the convective elements v. 
	Compare to sound speed. """

	delta_nabla_T = problem_1a()

	velocity_convective = (c.G * M / r**2 * delta_nabla_T / T)**(1/2) * l

	print "convective velocity is {0}".format(velocity_convective.decompose())

	velocity_sound = (gamma * P / rho)**(1/2)

	print "sound speed is {0}".format(velocity_sound.decompose())

	return (velocity_convective.decompose(), velocity_sound.decompose(), 
		(velocity_convective/velocity_sound).decompose())

def problem_1d():
	""" Reynolds number vl / nu """

	v = problem_1c()[0] 

	R = v * l / nu

	return R.decompose()

def problem_1e():
	""" The ratio of turbulent pressure to gas pressure."""

	v, v_s, x = problem_1c()

	return (v/v_s)**2

def problem_1f():
	""" The mean lifetime of the convective element t ~ l/v."""

	v = problem_1c()[0]

	convection_cell_lifetime = l / v

	print "Convection timescale: {0}".format(convection_cell_lifetime)

	kelvin_timescale = c.G * c.M_sun**2 / (c.R_sun * c.L_sun)

	print "Kelvin timescale: {0}".format(kelvin_timescale.decompose().to('yr'))

	nuclear_timescale = 0.007 * 0.3 * c.M_sun * c.c**2 / (c.L_sun)

	print "Nuclear timescale: {0}".format(nuclear_timescale.decompose().to('yr'))

	return convection_cell_lifetime.decompose()

