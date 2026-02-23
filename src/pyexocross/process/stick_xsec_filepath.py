def temperature_string_base(T, Tvib, Trot, NLTEMethod):
    """
    Base temperature string (no pressure/range), shared between stick spectra and cross sections.
    """
    if NLTEMethod == 'T' and Tvib is not None and Trot is not None:
        return f'Tvib{Tvib}K__Trot{Trot}K'
    elif NLTEMethod == 'D' and T is not None and Trot is not None:
        return f'T{T}K__Trot{Trot}K'
    else:
        # L or P (or fallback)
        T_val = T if T is not None else 0
        return f'T{T_val}K'


def temperature_pressure_string(T, P, temp_idx, NLTEMethod,
                                Tvib_list=None, Trot_list=None,
                                pressure_dependent=False):
    """
    Build temperature + pressure string for filenames.

    - Uses the same temperature form as stick spectra:
      * LTE / P:      T{T}K
      * Non-LTE T:    Tvib{Tvib}K__Trot{Trot}K
      * Non-LTE D:    T{T}K__Trot{Trot}K
    - Always uses the actual T / P for this file (no ranges).
    """
    # Temperature part
    if NLTEMethod == 'T' and temp_idx is not None and Tvib_list is not None and Trot_list is not None:
        Tvib_val = Tvib_list[temp_idx]
        Trot_val = Trot_list[temp_idx]
        T_str = temperature_string_base(T=None, Tvib=Tvib_val, Trot=Trot_val, NLTEMethod='T')
    else:
        T_str = temperature_string_base(T=T, Tvib=None, Trot=None, NLTEMethod=NLTEMethod)

    # Optional pressure part (single value, no range)
    if pressure_dependent and P is not None:
        if P < 0.001 or P >= 1000:
            P_str = f'__P{P:.2e}bar'
        else:
            P_str = f'__P{P}bar'
        T_str = T_str + P_str

    return T_str


def stick_spectra_filepath(ss_folder, T, Tvib, Trot, str_min_wnl, str_max_wnl, unit_fn,
                           data_info, wn_wl, UncFilter, threshold, database, abs_emi, LTE_NLTE, photo,
                           NLTEMethod):
    """
    Build stick spectra output file path (shared naming for ExoMol and HITRAN).
    """
    temp_part = temperature_string_base(T, Tvib, Trot, NLTEMethod)
    prefix = '__'.join(data_info) + '__' + temp_part + '__'
    return (ss_folder + prefix + wn_wl.lower() + str_min_wnl + '-' + str_max_wnl + unit_fn
            + 'unc' + str(UncFilter) + '__thres' + str(threshold) + '__' + database + '__'
            + abs_emi + photo + LTE_NLTE + '.stick')


def cross_section_filepath(xsecs_folder, data_info,
                           T, P, temp_idx,
                           Tvib_list, Trot_list,
                           str_min_v, str_max_v, unit_fn, wn_wl,
                           UncFilter, threshold, database, abs_emi,
                           bin_size, profile_label, LTE_NLTE, photo,
                           NLTEMethod, pressure_dependent):
    """
    Build cross-section (.xsec) output file path (shared naming for ExoMol and HITRAN).

    The temperature/pressure part follows the same rules as stick spectra and
    always uses the actual T / P for this file (no ranges).
    """
    temp_part = temperature_pressure_string(
        T=T,
        P=P,
        temp_idx=temp_idx,
        NLTEMethod=NLTEMethod,
        Tvib_list=Tvib_list,
        Trot_list=Trot_list,
        pressure_dependent=pressure_dependent,
    )
    prefix = '__'.join(data_info) + '__' + temp_part + '__'
    return (xsecs_folder + prefix + wn_wl.lower() + str_min_v + '-' + str_max_v + unit_fn
            + 'unc' + str(UncFilter) + '__thres' + str(threshold) 
            + '__BinSize' + str(bin_size) + unit_fn + profile_label.replace(' ', '')
            + '__' + database + '__' + abs_emi + photo + LTE_NLTE + '.xsec'
    )
