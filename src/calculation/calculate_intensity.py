"""
Calculate absorption and emission coefficients for LTE and non-LTE conditions.

This module provides functions for calculating line intensities (absorption/emission
coefficients) under various thermodynamic conditions.
"""
import numpy as np
import numexpr as ne
from ..base.constants import Inv8Pic, c2

## LTE intensity
# Calculate absorption coefficient
def cal_abscoefs(T_list, Q_list, Epp, gp, A, v, abundance):
    """
    Calculate the LTE absorption intensity for multiple temperatures.

    abscoef = gp * A * exp(- c2 * Epp / T) * (1 - exp(- c2 * v / T)) / (8 * pi * c * v^2 * Q) * abundance  

    Parameters
    ----------
    T_list : sequence of float
        Temperature list, shape (n_temps,)
    Q_list : sequence of float
        Partition function list, shape (n_temps,)
    Epp : np.ndarray
        Lower state energy array, shape (n_levels,)
    gp : np.ndarray
        Upper state degeneracy array, shape (n_levels,)
    A : np.ndarray
        Einstein A coefficient array, shape (n_levels,)
    v : np.ndarray
        Wavenumber array, shape (n_levels,)
    abundance : float
        Fractional abundance value
    
    Returns
    -------
    np.ndarray
        Absorption coefficient 2D array, shape (n_temps, n_levels)
    """
    # Ensure numpy arrays
    T_arr = np.asarray(T_list)[:, None]  
    Q_arr = np.asarray(Q_list)[:, None] 
    Epp_arr = np.asarray(Epp)[None, :]    
    gp_arr  = np.asarray(gp)[None, :]
    A_arr   = np.asarray(A)[None, :]
    v_arr   = np.asarray(v)[None, :]   
    absconst = ne.evaluate('gp_arr * A_arr * Inv8Pic / (v_arr ** 2) * abundance') 
    abscoef_arr = ne.evaluate('exp(- c2 * Epp_arr / T_arr) * (1 - exp(- c2 * v_arr / T_arr)) / Q_arr * absconst')   
    return abscoef_arr

# Calculate non-LTE absorption coefficient with two temperatures
def cal_abscoefs_nlte_2T(Tvib_list, Trot_list, Qnlte_arr, Evibpp, Erotpp, gp, A, v, abundance):
    """
    Calculate the non-LTE absorption intensity for multiple (Tvib, Trot) temperature pairs.

    abscoef_nlte =  ne.evaluate('gp * A * exp(- c2 * Evibpp / Tvib) * exp(- c2 * Erotpp / Trot) * (1 - exp(- c2 * v / Tvib)) / (8 * pi * c * v^2 * Q_nlte) * abundance')  
    

    Parameters
    ----------
    Tvib_list : sequence of float
        Vibrational temperature list, shape (n_temps,)
    Trot_list : sequence of float
        Rotational temperature list, shape (n_temps,)
    Qnlte_arr : np.ndarray
        Partition function array, shape (n_temps,)
    Evibpp : np.ndarray
        Lower state vibrational energy array, shape (n_levels,)
    Erotpp : np.ndarray
        Lower state rotational energy array, shape (n_levels,)
    gp : np.ndarray
        Upper state degeneracy array, shape (n_levels,)
    A : np.ndarray
        Einstein A coefficient array, shape (n_levels,)
    v : np.ndarray
        Wavenumber array, shape (n_levels,)
    abundance : float
        Fractional abundance value    
    
    Returns
    -------
    np.ndarray
        Non-LTE absorption coefficient 2D array, shape (n_temps, n_levels)
    """
    # Ensure numpy arrays
    Tvib_arr = np.asarray(Tvib_list)[:, None]  
    Trot_arr = np.asarray(Trot_list)[:, None]  
    Evibpp_arr = np.asarray(Evibpp)[None, :]
    Erotpp_arr = np.asarray(Erotpp)[None, :]
    gp_arr  = np.asarray(gp)[None, :]
    A_arr   = np.asarray(A)[None, :]
    v_arr   = np.asarray(v)[None, :]   
    absconst_nlte = ne.evaluate('gp_arr * A_arr * Inv8Pic / (v_arr ** 2) * abundance') 
    exponent = ne.evaluate('-c2 * (Evibpp_arr / Tvib_arr + Erotpp_arr / Trot_arr)')
    abscoef_nlte_arr = ne.evaluate('exp(exponent) * (1 - exp(- c2 * v_arr / Tvib_arr)) / Qnlte_arr * absconst_nlte')
    return abscoef_nlte_arr

# Calculate non-LTE absorption coefficient with density
def cal_abscoefs_nlte_nvib(T_list, Trot_list, Qnlte_arr, nvib, Erotpp, gp, A, v, abundance):
    """
    Calculate the non-LTE absorption intensity with vibrational density for multiple Trot.

    abscoef_nlte = ne.evaluate('gp * A * nvib * exp(- c2 * Erotpp / Trot) * (1 - exp(- c2 * v / T)) / (8 * pi * c * v^2 * Q_nlte) * abundance')

    Parameters
    ----------
    T_list : sequence of float
        Temperature list, shape (n_temps,)
    Trot_list : sequence of float
        Rotational temperature list, shape (n_temps,)
    Qnlte_arr : np.ndarray
        Partition function array, shape (n_temps,)
    nvib : np.ndarray
        Vibrational density array, shape (n_levels,)
    Erotpp : np.ndarray
        Lower state rotational energy array, shape (n_levels,)
    gp : np.ndarray
        Upper state degeneracy array, shape (n_levels,)
    A : np.ndarray
        Einstein A coefficient array, shape (n_levels,)
    v : np.ndarray
        Wavenumber array, shape (n_levels,)
    abundance : float
        Fractional abundance value    
    
    Returns
    -------
    np.ndarray
        Non-LTE absorption coefficient 2D array, shape (n_temps, n_levels)
    """
    # Ensure numpy arrays
    T_arr = np.asarray(T_list)[:, None]  
    Trot_arr = np.asarray(Trot_list)[:, None]  
    nvib_arr = np.asarray(nvib)[None, :]
    Erotpp_arr = np.asarray(Erotpp)[None, :]
    gp_arr  = np.asarray(gp)[None, :]
    A_arr   = np.asarray(A)[None, :]
    v_arr   = np.asarray(v)[None, :]   
    absconst_nlte = ne.evaluate('gp_arr * A_arr * Inv8Pic / (v_arr ** 2) * nvib_arr * abundance') 
    exponent = ne.evaluate('-c2 * (Erotpp_arr / Trot_arr)')
    abscoef_nlte_arr = ne.evaluate('exp(exponent) * (1 - exp(- c2 * v_arr / T_arr)) / Qnlte_arr * absconst_nlte')
    return abscoef_nlte_arr

### Using custom population
# Calculate non-LTE absorption coefficient with population
def cal_abscoefs_nlte_pop(T_list, pop, gp, gpp, A, v, abundance):
    """
    Calculate the non-LTE absorption intensity with non-LTE population for multiple T.

    abscoef_nlte = ne.evaluate('gp / gpp * A * pop * (1 - exp(- c2 * v / T)) / (8 * pi * c * v^2) * abundance')

    Parameters
    ----------
    T_list : sequence of float
        Temperature list, shape (n_temps,)
    pop : np.ndarray
        Non-LTE population array, shape (n_levels,)
    gp : np.ndarray
        Upper state degeneracy array, shape (n_levels,)
    gpp : np.ndarray
        Lower state degeneracy array, shape (n_levels,)
    A : np.ndarray
        Einstein A coefficient array, shape (n_levels,)
    v : np.ndarray
        Wavenumber array, shape (n_levels,)
    abundance : float
        Fractional abundance value    
    
    Returns
    -------
    np.ndarray
        Non-LTE absorption coefficient 2D array, shape (n_temps, n_levels)
    """
    # Ensure numpy arrays
    T_arr = np.asarray(T_list)[:, None]  
    pop_arr = np.asarray(pop)[None, :]
    gp_arr  = np.asarray(gp)[None, :]
    gpp_arr = np.asarray(gpp)[None, :]
    A_arr   = np.asarray(A)[None, :]
    v_arr   = np.asarray(v)[None, :]   
    absconst_nlte = ne.evaluate('gp_arr / gpp_arr * A_arr * Inv8Pic / (v_arr ** 2) * pop_arr * abundance') 
    abscoef_nlte_arr = ne.evaluate('(1 - exp(- c2 * v_arr / T_arr)) * absconst_nlte')
    return abscoef_nlte_arr
