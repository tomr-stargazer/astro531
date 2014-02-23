"""
The standard model is an approximation in which radiative, 
homogeneous stars are modeled by a polytrope of index n=3.
Using the standard model approximation, calculate the structure of the sun,
i.e. P(r) and rho(r). 
Chuck Cowley kindly provided the Lane-Emden function of index n=3 and 
you can find it in the file "Polytrope_n3.txt". 
Hint: Calculate the value of beta with at least 7 decimals; I found:
beta = 0.9995853.

The file bs05op.dat contains the "Standard Solar Model". 
This model is the preferred model for the interior of the sun 
(it does not include the atmosphere), 
calculated by Bahcall et al. 2005.
It reproduces the depth of the solar convection zone as inferred from 
heliseismological data and the measured neutrino flux.
The columns are self-explanatory.

Overplot the standard model values of P(r) and rho(r) on the values
of P(r) and rho(r) from the "Standard Solar Model".
Comment on the goodness of the polytropic approximation for the sun.

"""

from __future__ import division

import numpy as np

import scipy.optimize

import astropy
import astropy.units as u
import astropy.constants as c

from astro531.hw2.problem5 import problem_5a as hw2_problem_5a

# read in table
polytrope_n3 = astropy.table.Table.read("Polytrope_n3.txt", format='ascii',
                                        data_start=2)
polytrope_n3.rename_column('col1', 'Xi')
polytrope_n3.rename_column('col2', 'Theta')

def beta_solver(mass, mean_molecular_weight):

    mu = mean_molecular_weight

    # minimize this: 18 (1-beta)^1/2 / (mu^2 beta^2) - M = 0

    def beta_function(beta):
        return np.abs(18*( (1 - beta)**(1/2) / (mu**2 * beta**2) ) - mass)

    opt = scipy.optimize.minimize_scalar(beta_function, method='Bounded',
                                         bounds=[0,1], tol=1e-8)

    if opt.success:
        return opt.x
    else:
        raise Exception(opt.message)

def mu_from_abundances():
    """
    Mean molecular mass (including contribution from free electrons).

    Computed as 
    1/mu = X n_H/a_H + Y n_He/a_He + Z n_Z/a_Z

    where n_H/a_H is the number of particles per each unit mass of 
    hydrogen, and so on. Each hydrogen gives two particles, each
    helium (mass=4) gives three particles, and each metal (mass=N) 
    gives about N/2 particles.
    
    """

    XYZ = hw2_problem_5a()
    X = XYZ['X']
    Y = XYZ['Y']
    Z = 1 - X - Y

    one_over_mu = 2*X + 3/4*Y + 1/2*Z

    mu = 1/one_over_mu
    return mu

def polytrope_pressure_density(mean_molecular_weight=mu_from_abundances()):
    """
    Returns the pressure and density at all values of radius parameter "xi".

    Makes use of the polytrope_n3 table.

    """

    mu = mean_molecular_weight

    average_solar_density = c.M_sun / (4/3 * np.pi * c.R_sun**3)
    central_solar_density = 54.1825 * average_solar_density

    beta = beta_solver(1,1)

    density_array = central_solar_density * polytrope_n3['Theta']**3

    radiative_constant = 4 * c.sigma_sb / c.c
    a = radiative_constant

    k = ((c.k_B / (mu * c.u))**4 * 3/a * (1 - beta)/beta**4 )**(1/3)    

    pressure_array = k * density_array**(4/3)
    pressure_array_revised = (pressure_array.decompose().value * 
                              u.kg/(u.m * u.s**2))

    return pressure_array_revised.to('dyn cm-2'), density_array.to('g cm-3')    
    
    # from Cox, "Principles of Stellar Structure", Table 23.1
    #    central_pressure = 1.242e17 * u.dyn / (u.cm)**2
