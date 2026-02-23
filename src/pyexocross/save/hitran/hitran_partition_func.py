from pyexocross.save.exomol.exomol_partition_func import save_exomol_partition_func


def save_hitran_partition_func(states_df, Ntemp, Tmax):
    """
    Calculate and save partition functions for HITRAN database.

    Wrapper function that calls the ExoMol partition function calculation.

    Parameters
    ----------
    states_df : pd.DataFrame
        States DataFrame with 'E' and 'g' columns
    Ntemp : int
        Temperature step interval
    Tmax : int
        Maximum temperature in Kelvin
    """
    return save_exomol_partition_func(states_df, Ntemp, Tmax)
