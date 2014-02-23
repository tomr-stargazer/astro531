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

def polytrope_pressure_density(radius_parameter, mean_molecular_weight):
    """
    Returns the pressure and density at a given radius parameter "xi".

    Makes use of the polytrope_n3 table; throws errors when you try 
    to look up invalid values.

    """

    xi = radius_parameter
    mu = mean_molecular_weight

    try:
        theta = polytrope_n3['Theta'][polytrope_n3['Xi'] == xi][0]
    except IndexError:
        raise ValueError("Invalid input value of xi")

    average_solar_density = c.M_sun / (4/3 * np.pi * c.R_sun**3)

    central_solar_density = 54.1825 * average_solar_density

    # the second argument here SHOULD be the mean molecular weight of the Sun!
    beta = beta_solver(1,1)

    density_at_xi = central_solar_density * theta**3
    rho = density_at_xi

    radiative_constant = 4 * c.sigma_sb / c.c
    a = radiative_constant

    k = ((c.k_B / (mu * c.u))**4 * 3/a * (1 - beta)/beta**4 )**(1/3)
    #    k = 3.838e14 * (c.M_sun**(2/3) * c.R_sun**0)

    # from Cox, "Principles of Stellar Structure", Table 23.1
    #    central_pressure = 1.242e17 * u.dyn / (u.cm)**2

    pressure_at_xi = k * density_at_xi**(4/3)

    # unit hack because something isn't working right
    pressure_at_xi = pressure_at_xi.decompose().value * u.kg/(u.m * u.s**2)
    
    #    pressure_at_xi = central_pressure * theta**4
    
    return pressure_at_xi.to('dyn cm-2'), density_at_xi.to('g cm-3')
#    return pressure_at_xi.to('kg s-2 m'), density_at_xi

pressure = []
density = []

for xi in polytrope_n3['Xi']:

    P, rho = polytrope_pressure_density(xi, 0.7)
    pressure.append(P)
    density.append(rho)

    


    
