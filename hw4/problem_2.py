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

import astropy

# read in table
polytrope_n3 = astropy.table.Table.read("Polytrope_n3.txt", format='ascii',
                                        data_start=2)
polytrope_n3.rename_column('col1', 'Xi')
polytrope_n3.rename_column('col2', 'Theta')

def polytrope_pressure_density(radius_parameter):
    """
    Returns the pressure and density at a given radius parameter "xi".

    Makes use of the polytrope_n3 table; throws errors when you try 
    to look up invalid values.

    """

    xi = radius_parameter
    
