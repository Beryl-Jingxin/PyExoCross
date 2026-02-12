"""
Calculate emission coefficients for LTE and non-LTE conditions.

This module provides functions for calculating emission line intensities
under various thermodynamic conditions.
"""
import numpy as np
import numexpr as ne
from ..base.constants import hcInv4Pi, c2

## LTE emissivity
# Calculate emission coefficient
def cal_emicoefs(T_list, Q_list, Ep, gp, A, v, abundance):
    """
    Calculate the LTE emission intensity for multiple temperatures.

    emicoef = gp * h * c * A * v * exp(- c2 * Ep / T) / (4 * pi * Q) * abundance   

    Parameters
    ----------
    T_list : sequence of float
        Temperature list, shape (n_temps,)
    Q_list : sequence of float
        Partition function list, shape (n_temps,)
    Ep : np.ndarray
        Upper state energy array, shape (n_levels,)
    gp : np.ndarray
        Upper degeneracy array, shape (n_levels,)
    A : np.ndarray
        Einstein A coefficient array, shape (n_levels,)
    v : np.ndarray
        Wavenumber array, shape (n_levels,)
    abundance : float
        Fractional abundance value
    
    Returns
    -------
    np.ndarray
        Emission coefficient 2D array, shape (n_temps, n_levels)
    """
    # Ensure numpy arrays 
    T_arr = np.asarray(T_list)[:, None]  
    Q_arr = np.asarray(Q_list)[:, None] 
    Ep_arr = np.asarray(Ep)[None, :]  
    gp_arr  = np.asarray(gp)[None, :]
    A_arr   = np.asarray(A)[None, :]
    v_arr   = np.asarray(v)[None, :]    
    emiconst = ne.evaluate('gp_arr * A_arr * v_arr * hcInv4Pi * abundance') 
    emicoef_arr = ne.evaluate('exp(- c2 * Ep_arr / T_arr) / Q_arr * emiconst')   
    return emicoef_arr

# Calculate non-LTE emission coefficient with two temperatures
def cal_emicoefs_nlte_2T(Tvib_list, Trot_list, Qnlte_arr, Evibp, Erotp, gp, A, v, abundance):
    """
    Calculate the non-LTE emission intensity for multiple (Tvib, Trot) temperature pairs.

    emicoef_nlte = gp * A * v * h * c * exp(- c2 * Evibp / Tvib) * exp(- c2 * Erotp / Trot) / (4 * pi * Q_nlte) * abundance

    Parameters
    ----------
    Tvib_list : sequence of float
        Vibrational temperature list, shape (n_temps,)
    Trot_list : sequence of float
        Rotational temperature list, shape (n_temps,)
    Qnlte_arr : np.ndarray
        Partition function array, shape (n_temps,)
    Evibp : np.ndarray
        Upper state vibrational energy array, shape (n_levels,)
    Erotp : np.ndarray
        Upper state rotational energy array, shape (n_levels,)
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
        Non-LTE emission coefficient 2D array, shape (n_temps, n_levels)
    """
    # Ensure numpy arrays
    Tvib_arr = np.asarray(Tvib_list)[:, None]  
    Trot_arr = np.asarray(Trot_list)[:, None]  
    Evibp_arr = np.asarray(Evibp)[None, :]
    Erotp_arr = np.asarray(Erotp)[None, :]
    gp_arr  = np.asarray(gp)[None, :]
    A_arr   = np.asarray(A)[None, :]
    v_arr   = np.asarray(v)[None, :]   
    emiconst_nlte = ne.evaluate('gp_arr * A_arr * v_arr * hcInv4Pi * abundance') 
    exponent = ne.evaluate('-c2 * (Evibp_arr / Tvib_arr + Erotp_arr / Trot_arr)')
    emicoef_nlte_arr = ne.evaluate('exp(exponent) / Qnlte_arr * emiconst_nlte')
    return emicoef_nlte_arr

# Calculate non-LTE emission coefficient with density
def cal_emicoefs_nlte_nvib(Trot_list, Qnlte_arr, nvib, Erotp, gp, A, v, abundance):
    """
    Calculate the non-LTE emission intensity with vibrational density for multiple Trot.

    emicoef_nlte = gp * A * v * h * c * nvib * exp(- c2 * Erotp / Trot) / (4 * pi * Q_nlte) * abundance

    Parameters
    ----------
    Trot_list : sequence of float
        Rotational temperature list, shape (n_temps,)
    Qnlte_arr : np.ndarray
        Partition function array, shape (n_temps,)
    nvib : np.ndarray
        Vibrational density array, shape (n_levels,)
    Erotp : np.ndarray
        Upper state rotational energy array, shape (n_levels,)
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
        Non-LTE emission coefficient 2D array, shape (n_temps, n_levels)
    """
    # Ensure numpy arrays
    Trot_arr = np.asarray(Trot_list)[:, None]  
    nvib_arr = np.asarray(nvib)[None, :]
    Erotp_arr = np.asarray(Erotp)[None, :]
    gp_arr  = np.asarray(gp)[None, :]
    A_arr   = np.asarray(A)[None, :]
    v_arr   = np.asarray(v)[None, :]
    emiconst_nlte = ne.evaluate('gp_arr * A_arr * v_arr * hcInv4Pi * nvib_arr * abundance') 
    emicoef_nlte_arr = ne.evaluate('exp(-c2 * (Erotp_arr / Trot_arr)) / Qnlte_arr * emiconst_nlte')
    return emicoef_nlte_arr

# Calculate non-LTE emission coefficient with population
def cal_emicoefs_nlte_pop(pop, A, v, abundance):
    """
    Calculate the non-LTE emission intensity with non-LTE population.

    emicoef_nlte = A * v * h * c / (4 * pi) * pop * abundance

    Parameters
    ----------
    pop : np.ndarray
        Non-LTE population array, shape (n_levels,)
    A : np.ndarray
        Einstein A coefficient array, shape (n_levels,)
    v : np.ndarray
        Wavenumber array, shape (n_levels,)
    abundance : float
        Fractional abundance value
    
    Returns
    -------
    np.ndarray
        Non-LTE emission coefficient 2D array, shape (1, n_levels)
    """
    emicoef_nlte_arr = ne.evaluate('A * v * hcInv4Pi * pop * abundance')[None, :]
    return emicoef_nlte_arr
    