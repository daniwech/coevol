"""
Functions for the analysis of simulations can be put here.
"""

import copy
import pandas

def computeFitnessTraitDependence(ecoCoEvol, i, p, p_range):
    """

    :param ecoCoEvol:
    :param i:
    :param p:
    :param p_range:
    :return:
    """

    ecoCoEvol_ = copy.deepcopy(ecoCoEvol)

    dfResults = pandas.DataFrame(index=p_range, columns=["F"])
    for V_i_p in p_range:
        ecoCoEvol.V_[i*ecoCoEvol_.k+p] = V_i_p
        dN_t = ecoCoEvol.dN_dt()
        print dN_t
        dfResults.ix[V_i_p]['F'] = dN_t[i]

    return dfResults