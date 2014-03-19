"""
Problem 6.

The stellar structure equations are a system of four non-linear 
first-order differential equations which must be solved by numerical
integration. The solution of these four equations is complicated,
but for illustration purposes we can obtain a partial solution of 
the problem.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import trapz

import astropy.table
import astropy.constants as c
import astropy.units as u


solar_model = astropy.table.Table.read('../hw4/bs05op_revised.dat', 
                                       format='ascii.fixed_width_two_line', 
                                       data_start=1)

pressure_array = solar_model['P'] * u.dyn / u.cm**2
density_array = solar_model['Rho'] * u.g / u.cm**3
radius_array = (solar_model['R/Rsun'] * c.R_sun).to('cm')

epsilon_array = u.Quantity(3.038e-33, u.erg / u.g / u.s) * solar_model['X']**2 * density_array * solar_model['T']**(4.5)

luminosity_integrand = (epsilon_array * 4*np.pi* radius_array**2).to('erg cm-1 s-1')
mass_integrand = (density_array * 4*np.pi* radius_array**2).to('g cm-1')

calculated_luminosity_array = (np.zeros_like(solar_model['L/Lsun']) * luminosity_integrand.unit * u.cm).to('erg / s')
calculated_mass_array = (np.zeros_like(solar_model['M/Msun']) * mass_integrand.unit * u.cm).to('g')

for i in range(len(solar_model)):

	L_i = trapz(luminosity_integrand[:i], radius_array[:i])
	M_i = trapz(mass_integrand[:i], radius_array[:i])

	calculated_luminosity_array[i] = L_i
	calculated_mass_array[i] = M_i


def problem_6a():
	pass

def problem_6b():

	fig1 = plt.figure()
	plt.plot((radius_array/c.R_sun).decompose().value, 
		     (calculated_mass_array/c.M_sun).decompose().value)

	plt.xlabel("Radius r/R_sun")
	plt.ylabel("Mass M/M_sun")
	plt.title("Problem 6b. Tom Rice")

	fig2 = plt.figure()
	plt.plot((radius_array/c.R_sun).decompose().value, 
		     (calculated_luminosity_array/calculated_luminosity_array[-1]).decompose().value)

	plt.xlabel("Radius r/R_sun")
	plt.ylabel("Luminosity L/L_sun")
	plt.title("Problem 6b. Tom Rice")	

	plt.plot([0,1], [0.9,0.9], '--')

	return fig1, fig2

def problem_6c():

	fig1 = plt.figure()

	plt.plot(solar_model['R/Rsun'], solar_model['X'])
	plt.plot(solar_model['R/Rsun'], solar_model['Y(He4)'])

	plt.text(0.6, 0.77, "X: H mass fraction", color='b', fontsize=18)
	plt.text(0.6, 0.3, "Y: He mass fraction", color='g', fontsize=18)

	plt.xlabel("Radius r/R_sun")	
	plt.title("Problem 6c. Tom Rice")

	return fig1

def problem_6de():

	fig1 = plt.figure()

	plt.plot(solar_model['R/Rsun'], epsilon_array.value)

	plt.xlabel("Radius r/R_sun")	
	plt.ylabel("Epsilon")
	plt.title("Problem 6d. Tom Rice")	

	mu_array = 1 / (3* solar_model['X'] + 0.5 * (solar_model['Y(He4)'] + solar_model['He3']) + 1)

	fig2 = plt.figure()
	plt.plot(solar_model['R/Rsun'], mu_array)

	plt.xlabel("Radius r/R_sun")	
	plt.ylabel("Mean molecular weight")
	plt.title("Problem 6e. Tom Rice")	

	return fig1, fig2