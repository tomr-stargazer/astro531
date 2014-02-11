"""
Code to make some plots for HW3 in Nuria's class.

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import astropy.table

fermidirac_table = astropy.table.Table.read("FermiDiracIntegrals.txt", format='ascii', data_start=1)

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

    x_array = np.linspace(0, 120, 150)

    fig = plt.figure()
    
    for alpha in alpha_list:

        y_array = occupation_index(alpha, x_array)

        plt.plot(x_array, y_array, label=r"$\alpha$="+str(alpha))

    plt.show()        

    #    plt.semilogy()
    
    plt.xlabel(r"$E/kT$")
    plt.ylabel(r"$P(p)$")

    leg = plt.legend()
    leg.get_frame().set_alpha(0.5)

    return fig

        

        

    
