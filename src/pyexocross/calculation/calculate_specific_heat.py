"""
Calculate specific heat capacities for molecular systems.

This module provides functions for calculating specific heat at constant pressure
from partition function derivatives.
"""
import numpy as np
import numexpr as ne
from ..base.constants import c2, R

# Calculate Specific Heat
def cal_specific_heat(En, gn, T):
    """
    Calculate the specific heat capacity at constant pressure.

    Cp = R * (pfpp / pf - (pfp / pf)^2) + 2.5 * R

    Parameters
    ----------
    En : np.ndarray
        State energy array, shape (n_levels,)
    gn : np.ndarray
        State degeneracy array, shape (n_levels,)
    T : float
        Temperature in Kelvin

    Returns
    -------
    float
        Specific heat capacity at constant pressure (J/(K mol))
    """
    pf = ne.evaluate('sum(gn * exp(-c2 * En / T))')  
    pfp = ne.evaluate('sum(gn * exp(-c2 * En / T) * (c2 * En / T))')
    pfpp = ne.evaluate('sum(gn * exp(-c2 * En / T) * (c2 * En / T) ** 2)')
    cp = ne.evaluate('R * (pfpp / pf - (pfp / pf)**2) + 2.5 * R') 
    return cp