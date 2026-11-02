"""
Calculate partition functions for LTE and non-LTE conditions.

This module provides functions for calculating partition functions under
various thermodynamic conditions.
"""
import numpy as np
import numexpr as ne
from ..base.constants import c2


def partition_values(pf_df, T_list, state_loader):
    """Use tabulated integer Q(T) values and calculate all other values."""
    import warnings

    table = pf_df[['T', 'Q']].copy()
    table[['T', 'Q']] = table[['T', 'Q']].apply(
        lambda column: np.asarray(column, dtype=float)
    )
    table = table.dropna(subset=['T', 'Q'])
    if table.empty:
        raise ValueError('No valid partition function values in the partition function file.')

    temperatures = np.asarray(T_list, dtype=float)
    if np.any(temperatures <= 0):
        raise ValueError('Temperatures must be greater than 0 K.')

    values = dict(zip(table['T'], table['Q']))
    missing = list(dict.fromkeys(
        T for T in temperatures
        if not float(T).is_integer() or T not in values
    ))
    if missing:
        Tmin = float(table['T'].min())
        Tmax = float(table['T'].max())
        outside = [T for T in missing if T < Tmin or T > Tmax]
        if outside:
            requested = ', '.join(f'{T:g}' for T in outside)
            warnings.warn(
                f'Requested temperature(s) {requested} K are outside the partition '
                f'function file range {Tmin:g}-{Tmax:g} K (Tmax={Tmax:g} K). '
                'Q(T) will be calculated directly from the available states.',
                RuntimeWarning,
                stacklevel=2,
            )

        loaded = state_loader()
        chunks = [loaded] if isinstance(loaded, tuple) else loaded
        calculated = dict.fromkeys(missing, 0.0)
        for energies, degeneracies in chunks:
            energies = np.asarray(energies, dtype=float)
            degeneracies = np.asarray(degeneracies, dtype=float)
            for T in missing:
                calculated[T] += float(cal_partition_func(energies, degeneracies, T))
        values.update(calculated)

    return np.asarray([values[T] for T in temperatures], dtype=float)

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

