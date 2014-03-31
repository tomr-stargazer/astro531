"""
problem_4b.py

Astro 531 HW #6, problem 4

Check how well the slopes predictd in part (a) fit the data. 
Plot log L/Lsun vs. log M/Msun using data from the article
"M/L relation of intermediate-mass stars" (Malkov 2007, MNRAS, 382, 1073).
What criteria have been used to select the sample used in this article and 
why is this appropriate? At what mass is the break between the two slopes?

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import astropy.table

path = "/Users/tsrice/Dropbox/Grad School/Courses/Stars_astro531/"

malkov_table = astropy.table.Table.read(path+"Malkov_vizier_votable.vot")

fig = plt.figure()

# The data
plt.errorbar(malkov_table['Mass'], malkov_table['logL'], 
	xerr=malkov_table['e_Mass'], yerr=malkov_table['e_logL'], fmt='k.')

# The predictions
mass_array = np.logspace(-1, 2, 50)
upper_luminosity_array = (mass_array)**3
lower_luminosity_array = (mass_array)**(5.46)

plt.plot(mass_array, 0.85+np.log10(upper_luminosity_array), 'b--')
plt.plot(mass_array, np.log10(lower_luminosity_array), 'r--')

plt.semilogx()

plt.text(9, 3, r"Upper MS: $L \propto M^3$", color='b', fontsize=14)
plt.text(1, -1, r"Lower MS: $L \propto M^{5.46}$", color='r', fontsize=14)

plt.xlabel(r"Stellar mass ($M_\odot$)")
plt.ylabel(r"log ($L/L_\odot$)")

plt.title("Astro 531 HW6 problem 4b. Tom Rice")

plt.ylim(-3,6)

plt.show()
