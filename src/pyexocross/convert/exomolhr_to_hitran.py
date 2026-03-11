"""
Convert ExoMolHR format to HITRAN format.

ExoMolHR CSV files already contain precomputed line-list data (nu, unc, A, S,
E", g', g", J', J", QNs).  This module maps those columns directly into the
HITRAN .par fixed-width format.
"""
import numpy as np
import pandas as pd

from ..base import Timer, ensure_dir, print_conversion_info, print_file_info
from ..base.large_file import save_large_txt
from ..calculation.calculate_intensity import cal_abscoefs
from .exomol_to_hitran import error_code


def _ensure_conversion_globals():
    """Import legacy-style configuration globals needed for conversion."""
    from pyexocross.core import (
        ConversionMinFreq,
        ConversionMaxFreq,
        ConversionUnc,
        ConversionThreshold,
        GlobalQNLabel_list,
        GlobalQNFormat_list,
        LocalQNLabel_list,
        LocalQNFormat_list,
        read_path,
        save_path,
        data_info,
        abundance,
        species_main_id,
        species_sub_id,
    )
    globals().update(
        dict(
            ConversionMinFreq=ConversionMinFreq,
            ConversionMaxFreq=ConversionMaxFreq,
            ConversionUnc=ConversionUnc,
            ConversionThreshold=ConversionThreshold,
            GlobalQNLabel_list=GlobalQNLabel_list,
            GlobalQNFormat_list=GlobalQNFormat_list,
            LocalQNLabel_list=LocalQNLabel_list,
            LocalQNFormat_list=LocalQNFormat_list,
            read_path=read_path,
            save_path=save_path,
            data_info=data_info,
            abundance=abundance,
            species_main_id=species_main_id,
            species_sub_id=species_sub_id,
        )
    )


# ------------------------------------------------------------------ #
# Quantum-number formatting helpers
# ------------------------------------------------------------------ #

def _format_qn_column(series, fmt_str):
    """Format a single QN column according to its C-style format string.

    Parameters
    ----------
    series : pd.Series
        Column values (may be int, float, or str).
    fmt_str : str
        C-style format, e.g. ``'%12s'``, ``'%3d'``, ``'%5.1f'``.

    Returns
    -------
    pd.Series of str
    """
    py_fmt = fmt_str.replace('%', '{: >') + '}'
    if 'd' in py_fmt or 'f' in py_fmt:
        return pd.to_numeric(series, errors='coerce').map(py_fmt.format)
    else:
        return series.astype(str).str.strip().map(py_fmt.format)


def _format_qn_block(df, labels, formats, suffix):
    """Concatenate formatted QN columns into a single 15-char padded string.

    Parameters
    ----------
    df : pd.DataFrame
        ExoMolHR DataFrame (columns have ``'`` or ``"`` suffixes).
    labels : list of str
        QN label basenames (e.g. ``['+/-', 'e/f', 'ElecState', ...]``).
    formats : list of str
        Corresponding C-style format strings.
    suffix : str
        ``"'"`` for upper state, ``'"'`` for lower state.

    Returns
    -------
    pd.Series of str – each element is 15-char right-justified.
    """
    parts = []
    for label, fmt in zip(labels, formats):
        col_name = label + suffix
        if col_name in df.columns:
            parts.append(_format_qn_column(df[col_name], fmt))
        else:
            # Column not present – fill with blanks matching format width
            width = int(''.join(c for c in fmt if c.isdigit()).split('.')[0]) if any(c.isdigit() for c in fmt) else 1
            parts.append(pd.Series([' ' * width] * len(df), index=df.index))
    if not parts:
        return pd.Series([' ' * 15] * len(df), index=df.index)
    combined = pd.concat(parts, axis=1).astype(str).sum(axis=1)
    return combined.map('{: >15}'.format)


# ------------------------------------------------------------------ #
# Main conversion
# ------------------------------------------------------------------ #

def conversion_exomolhr2hitran(exomolhr_df):
    """Convert an ExoMolHR DataFrame to HITRAN ``.par`` format and save.

    Parameters
    ----------
    exomolhr_df : pd.DataFrame
        ExoMolHR line-list DataFrame as returned by ``read_exomolhr_df``.
    """
    from ..base.constants import Tref
    from ..database.load_exomolhr import read_exomolhr_pf

    _ensure_conversion_globals()

    print('Convert data format from ExoMolHR to HITRAN.')
    print_conversion_info(
        ConversionMinFreq, ConversionMaxFreq,
        GlobalQNLabel_list, GlobalQNFormat_list,
        LocalQNLabel_list, LocalQNFormat_list,
        ConversionUnc, ConversionThreshold,
    )

    t_tot = Timer()
    t_tot.start()

    # ---- Filter by frequency range ----
    df = exomolhr_df.copy()
    df = df[(df['v'] >= ConversionMinFreq) & (df['v'] <= ConversionMaxFreq)]

    # ---- Filter by uncertainty ----
    if ConversionUnc != 'None':
        df = df[df['unc'].astype(float) <= ConversionUnc]

    if len(df) == 0:
        raise ValueError(
            "No ExoMolHR transitions remain after applying conversion filters. "
            "Please check ConversionFrequncyRange and ConvUncFilter values."
        )

    # ---- Core arrays ----
    v = df['v'].values
    A = df['A'].values
    Epp = df['E"'].values
    gp = df["g'"].values.astype(int)
    gpp = df['g"'].values.astype(int)

    # ---- Intensity at Tref ----
    Q_ref = read_exomolhr_pf(read_path, data_info, [Tref])
    I = cal_abscoefs([Tref], Q_ref, Epp, gp, A, v, abundance)[0]

    # ---- Threshold filter ----
    if ConversionThreshold != 'None':
        mask = I >= ConversionThreshold
        df = df[mask].copy()
        v = v[mask]
        A = A[mask]
        Epp = Epp[mask]
        gp = gp[mask]
        gpp = gpp[mask]
        I = I[mask]

    if len(df) == 0:
        raise ValueError(
            "No ExoMolHR transitions remain after threshold filter. "
            "Please check ConvThreshold value."
        )

    nrows = len(df)

    # ---- Uncertainty → HITRAN error codes ----
    unc_values = df['unc'].values.astype(float).copy()
    unc_codes = error_code(unc_values)

    # ---- Default broadening (no .broad files for ExoMolHR) ----
    default_gamma = 0.07
    default_n_air = 0.50
    gamma_air = np.full(nrows, default_gamma)
    gamma_self = np.full(nrows, default_gamma)
    n_air = np.full(nrows, default_n_air)
    delta_air = [''] * nrows
    iref = [''] * nrows
    flag = [''] * nrows

    # ---- Quantum numbers ----
    Vp = _format_qn_block(df, GlobalQNLabel_list, GlobalQNFormat_list, "'")
    Vpp = _format_qn_block(df, GlobalQNLabel_list, GlobalQNFormat_list, '"')
    Qp = _format_qn_block(df, LocalQNLabel_list, LocalQNFormat_list, "'")
    Qpp = _format_qn_block(df, LocalQNLabel_list, LocalQNFormat_list, '"')

    QN_df = pd.DataFrame({
        "V'": Vp.values,
        'V"': Vpp.values,
        "Q'": Qp.values,
        'Q"': Qpp.values,
    })

    # ---- Assemble HITRAN DataFrame ----
    hitran_df = pd.DataFrame({
        'M': str(species_main_id),
        'I': str(species_sub_id),
        'v': v,
        'S': I,
        'A': A,
        'gamma_air': gamma_air,
        'gamma_self': gamma_self,
        'E"': Epp,
        'n_air': n_air,
        'delta_air': delta_air,
    })
    hitran_df = pd.concat([hitran_df, QN_df], axis=1)
    hitran_df['Ierr'] = unc_codes
    hitran_df['Iref'] = iref
    hitran_df['flag'] = flag
    hitran_df["g'"] = gp.astype(float)
    hitran_df['g"'] = gpp.astype(float)

    # Sort by wavenumber
    hitran_df.sort_values('v', inplace=True)

    t_tot.end()
    print('Finished converting ExoMolHR linelist to HITRAN format!\n')

    # ---- Save ----
    print('Saving HITRAN format data into file ...')
    ts = Timer()
    ts.start()
    conversion_folder = save_path + 'conversion/ExoMolHR2HITRAN/'
    ensure_dir(conversion_folder)
    conversion_path = conversion_folder + '__'.join(data_info[-2:]) + '.par'
    hitran_format = (
        "%2s%1s%12.6f%10.3E%10.3E%5.3f%5.3f%10.4f%4.2f%8s"
        "%15s%15s%15s%15s%6s%12s%1s%7.1f%7.1f"
    )
    save_large_txt(conversion_path, hitran_df, fmt=hitran_format)
    ts.end()

    hitran_fmt_list = ['%' + i for i in hitran_format.split('%')][1:]
    print_file_info('Converted HITRAN par', hitran_df.columns, hitran_fmt_list)
    print('Converted HITRAN par file has been saved:', conversion_path)
    print('Converted HITRAN par file has been saved!\n')
    print('Finished converting data format from ExoMolHR to HITRAN!\n')
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')
