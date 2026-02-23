"""
Calculate oscillator strengths from Einstein A coefficients.

This module provides functions for converting between Einstein A coefficients
and oscillator strengths (gf or f values).
"""
import numpy as np
import numexpr as ne
from ..base.constants import c

# Calculate Oscillator Strength
# Calculate gf or f
def cal_oscillator_strength(gp, gpp, A, v, gfORf='f'):
    """
    Calculate oscillator strength (gf or f) from Einstein A coefficient.

    gf = gp * A / (c * v)^2
    f = gf / gpp

    Parameters
    ----------
    gp : np.ndarray
        Upper state degeneracy array, shape (n_levels,)
    gpp : np.ndarray
        Lower state degeneracy array, shape (n_levels,)
    A : np.ndarray
        Einstein A coefficient array, shape (n_levels,)
    v : np.ndarray
        Wavenumber array, shape (n_levels,)

    Returns
    -------
    np.ndarray
        Oscillator strength array (gf or f depending on global gfORf setting),
        shape (n_levels,)
    """
    gf = ne.evaluate('gp * A / (c * v)**2')
    if 'G' in gfORf.upper():
        return gf
    else:
        f = ne.evaluate('gf / gpp')
        return f
