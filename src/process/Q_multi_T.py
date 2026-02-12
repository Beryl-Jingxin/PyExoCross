"""
Calculate partition functions for multiple temperatures.

This module provides functions for calculating partition functions under
LTE and various non-LTE conditions.
"""
import numpy as np
import pandas as pd
from ..calculation.calculate_partition_func import cal_Q_nlte_2T, cal_Q_nlte_nvib
from ..database.load_exomol import read_exomol_pf

def cal_pf_multiT(T_list, Tvib_list, Trot_list, states_df, NLTEMethod, read_path, data_info):
    """
    Calculate partition functions for multiple temperatures using specified method.

    Supports LTE and various non-LTE methods (two-temperature, custom density,
    custom population).

    Parameters
    ----------
    T_list : list of float
        Temperature list for LTE method
    Tvib_list : list of float
        Vibrational temperature list for non-LTE two-temperature method
    Trot_list : list of float
        Rotational temperature list for non-LTE methods
    states_df : pd.DataFrame
        States DataFrame with energy and degeneracy columns, and optionally
        Evib, Erot, nvib columns depending on method

    Returns
    -------
    np.ndarray
        Partition function array, shape (n_temps,)
    """
    # LTE
    if NLTEMethod == 'L':
        Q_arr = read_exomol_pf(read_path, data_info, T_list)
    # Non-LTE using two temperatures
    elif NLTEMethod == 'T':
        Q_arr = cal_Q_nlte_2T(Tvib_list, Trot_list, states_df['Evib'], states_df['Erot'], states_df['g'])
    # Non-LTE using custom density
    elif NLTEMethod == 'D':
        Q_arr = cal_Q_nlte_nvib(Trot_list, states_df['nvib'], states_df['E'], states_df['g'])
    # Non-LTE using custom population   
    elif NLTEMethod == 'P':
        # Q is not used in population-based calculations, return array of ones with matching length
        Q_arr = np.ones(len(T_list), dtype=float)
    else:
        raise ValueError("Please choose one non-LTE method from: 'T', 'D' or 'P'.")    
    return(Q_arr)
