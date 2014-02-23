"""
Code to compute abundance values X, Y, Z (etc) for Problem 5 of 
Astro 531 homework #2.

"""

from __future__ import division

import numpy as np

import astropy.table

path = "/Users/tsrice/Documents/Code/astro531/hw2/"

cowley_table = astropy.table.Table.read(path+"Table3.2.Cowley.txt", format='ascii')

cowley_table.rename_column('col1', 'element')
cowley_table.rename_column('col2', 'symbol')
cowley_table.rename_column('col3', 'atomic_number')
cowley_table.rename_column('col4', 'atomic_weight')
cowley_table.rename_column('col5', 'log_abundance')

def un_dexed_abundance(dex_abundance):
    """ Gets the fractional abundance of an element from its dex. """
    return 10**(dex_abundance - 12)

total_mass = np.sum(
    un_dexed_abundance(cowley_table['log_abundance'])*cowley_table['atomic_weight'])

#print total_mass

# X: mass fraction of H.   
# Y: mass fraction of He.
# Z: mass fraction of "metals"

def problem_5a():

    X = (un_dexed_abundance(cowley_table['log_abundance'][0]) * 
         cowley_table['atomic_weight'][0]) / total_mass

    Y = (un_dexed_abundance(cowley_table['log_abundance'][1]) *
         cowley_table['atomic_weight'][1]) / total_mass

    Z = 1 - X - Y

    #    print "X: " + str(X)
    #    print "Y: " + str(Y)
    #    print "Z: " + str(Z)

    return {'X':X, 'Y':Y, 'Z':Z}

def calculate_abundance_relative_to_H(table, symbol):

    # What row of the table has this symbol's data?
    index, = np.where(table['symbol'] == symbol)[0]

    number_abundance = un_dexed_abundance(table['log_abundance'][index])
    mass_abundance = (number_abundance * table['atomic_weight'][index] /
                      table['atomic_weight'][0])

    return number_abundance, mass_abundance

assert calculate_abundance_relative_to_H(cowley_table, 'H') == (1., 1.)

def problem_5b():

    symbol_list_ii = ['C', 'N', 'O', 'Ne']
    symbol_list_iii = ['Cr', 'Mn', 'Fe', 'Co', 'Ni']

    print "element : number abundance : mass abundance"

    for symbol_list in [['He'], symbol_list_ii, symbol_list_iii]:
        for symbol in symbol_list:

            number_abundance, mass_abundance = \
              calculate_abundance_relative_to_H(cowley_table, symbol)
            print symbol," : ", number_abundance, " : ", mass_abundance

            
