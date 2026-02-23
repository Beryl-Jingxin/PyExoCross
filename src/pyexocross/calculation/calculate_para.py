"""
Calculate various molecular parameters from energy levels and transitions.

This module provides functions for calculating wavenumbers, energies, quantum numbers,
and other derived parameters from state and transition data.
"""
import numpy as np
import numexpr as ne
import pandas as pd
from ..base.constants import c2

# Convert frequency, upper and lower energy and J
# Calculate frequency
def cal_v(Ep, Epp):
    """
    Calculate transition wavenumber from upper and lower state energies.

    v = Ep - Epp

    Parameters
    ----------
    Ep : np.ndarray or float
        Upper state energy values
    Epp : np.ndarray or float
        Lower state energy values

    Returns
    -------
    np.ndarray or float
        Transition wavenumber(s) computed as Ep - Epp
    """
    v = ne.evaluate('Ep - Epp')
    return v

# Calculate upper state energy with ExoMol database
def cal_Ep(Epp, v):
    """
    Calculate upper state energy from lower state energy and wavenumber.

    Ep = Epp + v

    Parameters
    ----------
    Epp : np.ndarray or float
        Lower state energy values
    v : np.ndarray or float
        Transition wavenumber values

    Returns
    -------
    np.ndarray or float
        Upper state energy values
    """
    Ep = ne.evaluate('Epp + v')
    return Ep

# Calculate upper state energy with HITRAN database
def cal_Ep_hitran(hitran_df):
    """
    Calculate upper state energy from HITRAN DataFrame.

    Parameters
    ----------
    hitran_df : pd.DataFrame
        DataFrame containing 'Epp' and 'v' columns

    Returns
    -------
    pd.DataFrame
        DataFrame with 'Ep' column containing upper state energies
    """
    Epp = hitran_df['Epp'].values
    v = hitran_df['v'].values
    Ep = cal_Ep(Epp, v)
    Ep_df = pd.DataFrame(Ep,columns=['Ep'])
    return Ep_df

# Calculate upper J
def cal_Jp(Fp, Fpp, Jpp):
    """
    Calculate upper state total angular momentum quantum number.

    Jp = Fp + Fpp - Jpp

    Parameters
    ----------
    Fp : np.ndarray or float
        Upper state F quantum number
    Fpp : np.ndarray or float
        Lower state F quantum number
    Jpp : np.ndarray or float
        Lower state J quantum number

    Returns
    -------
    np.ndarray or float
        Upper state J quantum number
    """
    Jp = ne.evaluate('Fp + Fpp - Jpp')
    return Jp

# Calculate F
def cal_F(g):
    """
    Calculate F quantum number from degeneracy.

    F = (g - 1) / 2

    Parameters
    ----------
    g : np.ndarray or float
        Degeneracy values

    Returns
    -------
    np.ndarray or float
        F quantum number values
    """
    F = ne.evaluate('(g - 1) * 0.5')
    return F

# Calculate uncertainty
def cal_uncertainty(unc_u, unc_l):
    """
    Calculate combined uncertainty from upper and lower state uncertainties.

    unc = sqrt(unc_u^2 + unc_l^2)

    Parameters
    ----------
    unc_u : np.ndarray or float
        Upper state uncertainty values
    unc_l : np.ndarray or float
        Lower state uncertainty values

    Returns
    -------
    np.ndarray or float
        Combined uncertainty values
    """
    unc = ne.evaluate('sqrt(unc_u ** 2 + unc_l ** 2)')
    return unc