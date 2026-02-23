def get_ntemp(NLTEMethod, T_list, Trot_list):
    """
    Number of temperature indices for stick spectra / cross section loops.

    LTE: ntemp = len(T_list).
    Non-LTE T or D: ntemp = len(Trot_list).
    Non-LTE P: one output per temperature in T_list.

    Parameters
    ----------
    NLTEMethod : str
        'L', 'T', 'D', or 'P'
    T_list : list
        Temperature list (LTE)
    Trot_list : list
        Rotational temperature list (Non-LTE T/D)

    Returns
    -------
    int
        Number of iterations (temperature indices).
    """
    from pyexocross.core import abs_emi
    if NLTEMethod == 'L':
        return len(T_list)
    if NLTEMethod in ('T', 'D'):
        return len(Trot_list)
    if NLTEMethod == 'P':
        # For population-based non-LTE, users still expect one output per input
        # temperature for both absorption and emission calculations.
        return len(T_list)
    return 0


def get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list):
    """
    Get (T, Tvib, Trot) for a given temperature index.

    Parameters
    ----------
    temp_idx : int
        Temperature index (0 .. ntemp-1).
    NLTEMethod : str
        'L', 'T', 'D', or 'P'
    T_list : list
        Temperature list
    Tvib_list : list
        Vibrational temperature list (Non-LTE T only)
    Trot_list : list
        Rotational temperature list (Non-LTE T/D)

    Returns
    -------
    tuple of (float or None, float or None, float or None)
        (T, Tvib, Trot). Unused entries are None.
    """
    from pyexocross.core import abs_emi
    if NLTEMethod == 'L':
        T = T_list[temp_idx] if temp_idx < len(T_list) else None
        return (T, None, None)
    if NLTEMethod == 'T':
        Tvib = Tvib_list[temp_idx] if temp_idx < len(Tvib_list) else None
        Trot = Trot_list[temp_idx] if temp_idx < len(Trot_list) else None
        T = Trot
        return (T, Tvib, Trot)
    if NLTEMethod == 'D':
        T = T_list[temp_idx] if temp_idx < len(T_list) else None
        Trot = Trot_list[temp_idx] if temp_idx < len(Trot_list) else None
        return (T, None, Trot)
    if NLTEMethod == 'P':
        T = T_list[temp_idx] if temp_idx < len(T_list) else (T_list[0] if T_list else None)
        return (T, None, None)
    return (None, None, None)