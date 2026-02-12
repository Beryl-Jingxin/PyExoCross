from src.save.exomol.exomol_lifetime import save_exomol_lifetime


def save_hitran_lifetime(read_hitran2exomol_path, states_df, hitran_states_col, hitran_states_fmt):
    """
    Main function to calculate and save radiative lifetimes for HITRAN database.

    Wrapper function that calls the ExoMol lifetime calculation after converting
    HITRAN data to ExoMol format.

    Parameters
    ----------
    read_hitran2exomol_path : str
        Path to converted HITRAN-to-ExoMol transition files
    states_df : pd.DataFrame
        States DataFrame
    hitran_states_col : list of str
        Column names for HITRAN states file
    hitran_states_fmt : list of str
        Format strings for each column
    """
    return save_exomol_lifetime(read_hitran2exomol_path, states_df, hitran_states_col, hitran_states_fmt)
