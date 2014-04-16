"""
Problem 1 from Nuria's homework.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
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

t_array = np.logspace(11, 18, 100)
t_array_years = t_array / (31557600.0)

hr_array = hayashi_radius(t_array)
Tc_array = central_temperature(hr_array)

def plot_problem1d():
	fig1 = plt.figure()
	plt.plot(t_array_years, (hr_array/c.R_sun).value, lw=2)
	plt.title("Hayashi radius")
	plt.xlabel("Time (years)")
	plt.ylabel(r"Hayashi radius ($R_\odot$)")
	plt.xlim(t_array_years[0], t_array_years[-1])
	plt.semilogx()

	fig2 = plt.figure()
	plt.plot(t_array_years, Tc_array.value, lw=2)
	plt.title("Central temperature")
	plt.ylabel("Central temperature (K)")
	plt.xlabel("Time (years)")
	plt.xlim(t_array_years[0], t_array_years[-1])
	plt.semilogx()

	plt.show()

	return fig1, fig2


plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.minor.size'] = 1.5
plt.rcParams['ytick.minor.size'] = 1.5
plt.rcParams['xtick.major.size'] = 2.0
plt.rcParams['ytick.major.size'] = 2.0
font = {'family' : 'normal',
	'weight' : 'normal',
	'size' : 14}
matplotlib.rc('font', **font)

fig1, fig2 = plot_problem1d()