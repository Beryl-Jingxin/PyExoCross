"""
Calculate partition functions for LTE and non-LTE conditions.

This module provides functions for calculating partition functions under
various thermodynamic conditions.
"""
import numpy as np
import numexpr as ne
from ..base.constants import c2

## LTE partition function
# Calculate Parition Function
def cal_partition_func(En, gn, T):
    """
    Calculate the LTE partition function for a single temperature.

    Q(T) = sum(gn * exp(-c2 * En / T))

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
        Partition function value
    """
    pf = ne.evaluate('sum(gn * exp(-c2 * En / T))') 
    return pf
    
## Non-LTE partition function
### Using two temperatures
# Calculate non-LTE partition function
def cal_Q_nlte_2T(Tvib_list, Trot_list, Evib, Erot, g):
    """
    Calculate the non-LTE partition function for multiple (Tvib, Trot) temperature pairs.

    Q_nlte = sum(g * exp(-c2 * Evib / Tvib) * exp(-c2 * Erot / Trot))

    Parameters
    ----------
    Tvib_list : sequence of float
        Vibrational temperature list, shape (n_temps,)
    Trot_list : sequence of float
        Rotational temperature list, shape (n_temps,)
    Evib : np.ndarray
        State vibrational energy array, shape (n_levels,)
    Erot : np.ndarray
        State rotational energy array, shape (n_levels,)
    g : np.ndarray
        State degeneracy array, shape (n_levels,)  
    
    Returns
    -------
    np.ndarray
        Non-LTE partition function array, shape (n_temps,)
    """
    Tvib_arr = np.asarray(Tvib_list)[:, None]   
    Trot_arr = np.asarray(Trot_list)[:, None]       
    Evib_arr = np.asarray(Evib)[None, :]   
    Erot_arr = np.asarray(Erot)[None, :]   
    g_arr = np.asarray(g)[None, :]  
    exponent = ne.evaluate('-c2 * (Evib_arr / Tvib_arr + Erot_arr / Trot_arr)')
    Qnlte_arr = (g_arr * np.exp(exponent)).sum(axis=1)
    return Qnlte_arr

### Using custom density
# Calculate non-LTE partition function
def cal_Q_nlte_nvib(Trot_list, nvib, Erot, g):
    """
    Calculate the non-LTE partition function for multiple Trot.

    Q_nlte = sum(g * nvib * exp(-c2 * Erot / Trot))
    
    Parameters
    ----------
    Trot_list : sequence of float
        Rotational temperature list, shape (n_temps,)
    nvib : np.ndarray
        Vibrational population array, shape (n_levels,)
    Erot : np.ndarray
        Rotational energy array, shape (n_levels,)
    g : np.ndarray
        Degeneracy array, shape (n_levels,)
    
    Returns
    -------
    np.ndarray
        Non-LTE partition function array, shape (n_temps,)
    """
    # Ensure numpy arrays 
    Trot_arr = np.asarray(Trot_list)[:, None]    
    nvib_arr = np.asarray(nvib)[None, :]    
    Erot_arr = np.asarray(Erot)[None, :] 
    g_arr = np.asarray(g)[None, :]     
    exponent = ne.evaluate('-c2 * (Erot_arr / Trot_arr)')
    Qnlte_arr = (g_arr * nvib_arr * np.exp(exponent)).sum(axis=1)   
    return Qnlte_arr
    