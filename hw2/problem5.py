"""
Code to compute abundance values X, Y, Z (etc) for Problem 5 of 
Astro 531 homework #2.

"""

from __future__ import division

import numpy as np

import astropy

cowley_table = astropy.table.Table.read("Table3.2.Cowley.txt", format='ascii')

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

print total_mass

# X: mass fraction of H.   
# Y: mass fraction of He.
# Z: mass fraction of "metals"

X = (un_dexed_abundance(cowley_table['log_abundance'][0]) * 
     cowley_table['atomic_weight'][0]) / total_mass

Y = (un_dexed_abundance(cowley_table['log_abundance'][1]) *
     cowley_table['atomic_weight'][1]) / total_mass

Z = 1 - X - Y

print "X: " + str(X)
print "Y: " + str(Y)
print "Z: " + str(Z)
