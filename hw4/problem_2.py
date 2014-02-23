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
import matplotlib.pyplot as plt

import scipy.optimize

import astropy
import astropy.units as u
import astropy.constants as c

from astro531.hw2.problem5 import problem_5a as hw2_problem_5a

# read in polytrope table
polytrope_n3 = astropy.table.Table.read("Polytrope_n3.txt", format='ascii',
                                        data_start=2)
polytrope_n3.rename_column('col1', 'Xi')
polytrope_n3.rename_column('col2', 'Theta')

# read in solar table
solar_model = astropy.table.Table.read('bs05op_revised.dat', 
                                       format='ascii.fixed_width_two_line', 
                                       data_start=1)

def beta_solver(mass, mean_molecular_weight):
    """
    Solves for the beta radiation parameter.
    
    Beta is the parameter describing what proportion of the pressure
    is due to radiation versus gas pressure.

    """

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

    XYZ_dict = hw2_problem_5a()
    X = XYZ_dict['X']
    Y = XYZ_dict['Y']
    Z = 1 - X - Y

    one_over_mu = 2*X + 3/4*Y + 1/2*Z

    mu = 1/one_over_mu
    return mu

def polytrope_pressure_density_radius(mean_molecular_weight=mu_from_abundances()):
    """
    Returns the pressure and density at all values of radius parameter "xi".

    Makes use of the polytrope_n3 table.

    """

    mu = mean_molecular_weight

    # bulk density equals mass over volume
    average_solar_density = c.M_sun / (4/3 * np.pi * c.R_sun**3)
    # from Cox, "Principles of Stellar Structure", Table 23.1    
    central_solar_density = 54.1825 * average_solar_density

    beta = beta_solver(1,mu)

    # from Cox, "Principles of Stellar Structure", Table 23.1    
    xi_at_stellar_surface = 6.89685
    radius_array = c.R_sun * polytrope_n3['Xi'] / xi_at_stellar_surface
    
    density_array = central_solar_density * polytrope_n3['Theta']**3

    # radiative_constant 
    a_rad = 4 * c.sigma_sb / c.c

    k = ((c.k_B / (mu * c.u))**4 * 3/a_rad * (1 - beta)/beta**4 )**(1/3)

    pressure_array = k * density_array**(4/3)

    print pressure_array[0]
    pressure_array_revised = (pressure_array.decompose().value * 
                              u.kg/(u.m * u.s**2))

    return pressure_array_revised.to('dyn cm-2'), density_array.to('g cm-3'), radius_array.to('cm')
    
    # from Cox, "Principles of Stellar Structure", Table 23.1
    #    central_pressure = 1.242e17 * u.dyn / (u.cm)**2

def solar_model_pressure_density_radius():
    """
    Extracts the solar pressure & density from the attached table.

    """

    pressure_array = solar_model['P'] * u.dyn / u.cm**2
    density_array = solar_model['Rho'] * u.g / u.cm**3
    radius_array = solar_model['R/Rsun'] * c.R_sun

    return pressure_array, density_array, radius_array.to('cm')
    
    
def plot_polytrope_and_solar_model():

    poly_P, poly_rho, poly_r = polytrope_pressure_density_radius()
    solar_P, solar_rho, solar_r = solar_model_pressure_density_radius()

    # Plot P(r) 
    fig1 = plt.figure()

    plt.plot(poly_r, poly_P, label='n=3 polytrope')
    plt.plot(solar_r, solar_P, label='bs05 standard solar model')

    plt.legend()

    plt.xlabel("Radius (cm)")
    plt.ylabel(r"Pressure (dyn/cm$^2$)")

    plt.title("Pressure vs radius in polytrope and standard solar model.")
               
    # Plot rho(r)
    fig2 = plt.figure()

    plt.plot(poly_r, poly_rho, label='n=3 polytrope')
    plt.plot(solar_r, solar_rho, label='bs05 standard solar model')

    plt.legend()

    plt.xlabel("Radius (cm)")
    plt.ylabel(r"density (g/cm$^3$)")

    plt.title("Density vs radius in polytrope and standard solar model.")    
    
    return fig1, fig2
