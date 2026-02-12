from src.save.exomol.exomol_specific_heat import save_exomol_specific_heat


def save_hitran_specific_heat(states_df, Ntemp, Tmax):
    """
    Calculate and save specific heat capacities for HITRAN database.

    Wrapper function that calls the ExoMol specific heat calculation.

    Parameters
    ----------
    states_df : pd.DataFrame
        States DataFrame with 'E' and 'g' columns
    Ntemp : int
        Temperature step interval
    Tmax : int
        Maximum temperature in Kelvin
    """
    return save_exomol_specific_heat(states_df, Ntemp, Tmax)
    