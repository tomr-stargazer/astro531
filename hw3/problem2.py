"""
Code to make some plots for HW3 in Nuria's class.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import astropy.table
import astropy.constants as c
import astropy.units as u
from astropy.units.quantity import Quantity

fermidirac_table = astropy.table.Table.read("FermiDiracIntegrals.txt", 
                                            format='ascii', data_start=1)

fermidirac_table.rename_column('col1', 'alpha')
fermidirac_table.rename_column('col2', '2/3 F_3/2')
fermidirac_table.rename_column('col3', 'F_1/2')


def occupation_index(alpha, energy):
    """ Calculates the occupation index P(p). """

    P_p = 1/(np.exp(alpha + energy) + 1)

    return P_p

def problem_2b():
    """
    Calculate and plot the occupation index P(p) 
    as a function of E/kT for alpha = 10, 0, -4, -30, -100.

    """

    alpha_list = [10, 0, -4, -30, -100]

    x_array = np.linspace(-20, 120, 150)

    fig = plt.figure()
    
    for alpha in alpha_list:

        y_array = occupation_index(alpha, x_array)

        plt.plot(x_array, y_array, label=r"$\alpha$="+str(alpha))

    plt.show()        

    plt.xlabel(r"$E/kT$")
    plt.ylabel(r"$P(p)$")

    leg = plt.legend()
    leg.get_frame().set_alpha(0.5)

    return fig

def electron_number_density(alpha, temperature):
    """
    Calculates the electron number density given an input alpha & T.

    """

    T = Quantity(temperature, u.K)

    F_12_alpha = fermidirac_table['F_1/2'][alpha == fermidirac_table['alpha']]

    n_e = 4 * np.pi / c.h**3 * (2 * c.m_e * c.k_B * T)**(3/2) * F_12_alpha

    return n_e.to('1/cm3')

def electron_pressure(alpha, temperature):

    T = Quantity(temperature, u.K)

    n_e = electron_number_density(alpha, temperature)

    twothirds_F_32_alpha = fermidirac_table['2/3 F_3/2'][alpha == 
                                                  fermidirac_table['alpha']]
    F_12_alpha = fermidirac_table['F_1/2'][alpha == fermidirac_table['alpha']]    
    P_e = n_e * c.k_B * T * (twothirds_F_32_alpha / F_12_alpha)

    return P_e.to('dyn / cm2')
  

def problem_2c():
    """
    The attached table gives the Fermi-Dirac integrals for a range of 
    values of the parameter alpha.
    Assuming that the pressure is mostly due to electrons, calculate 
    the pressure P, electron density, and density rho for all values 
    of alpha in the table and temperatures 10^7 K, 10^8 K, and 10^9 K. 
    Plot P vs n_e for each temperature.
    Indicate the value of alpha along the curve.
    
    """

    fig = plt.figure()

    alpha_array = fermidirac_table['alpha']
    T_array = np.array([1e9, 1e8, 1e7]) * u.K

    for T in T_array:

        n_e_array = electron_number_density(alpha_array, T)
        
        P_e_array = electron_pressure(alpha_array, T)

        plt.plot(n_e_array, P_e_array, label=r"T = 10$^%s$" % 
                 str(int(np.log10(T.value)))+" K")

    plt.xlabel(r"$n_e$ (cm$^{-2}$)")
    plt.ylabel(r"$P_e$ (dyn cm$^{-2}$)")

    plt.loglog()

    plt.legend(loc='upper left')

    return fig

        
    
