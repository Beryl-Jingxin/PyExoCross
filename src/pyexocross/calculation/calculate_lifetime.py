"""
Calculate radiative lifetimes from transition probabilities.

This module provides functions for calculating radiative lifetimes of excited states
from Einstein A coefficients.
"""
import numpy as np
import pandas as pd

# Calculate Lifetime
def cal_lifetime(states_df, trans_df):
    """
    Calculate radiative lifetime for each state from transition probabilities.

    Lifetime is computed as the inverse of the sum of Einstein A coefficients
    for all transitions originating from each upper state.

    Parameters
    ----------
    states_df : pd.DataFrame
        DataFrame containing state information with 'id' column.
    trans_df : pd.DataFrame
        DataFrame containing transitions with 'u' (upper state) and 'A' (Einstein A) columns.

    Returns
    -------
    np.ndarray
        Lifetime array corresponding to states in states_df, shape (n_states,)
    """
    A_sum = trans_df.groupby('u')['A'].sum()
    ids = states_df['id']
    lifetime_whole = np.zeros(max(ids)+1)
    lifetime_whole[A_sum.index] = A_sum.values
    lt = lifetime_whole[ids]
    return lt
    