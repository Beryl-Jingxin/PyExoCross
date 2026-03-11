"""
Load ExoMolHR database files.

This module provides helpers for reading ExoMolHR CSV line lists and
partition-function files using a line-list workflow similar to HITRAN.
"""
import glob
import os

import numpy as np
import pandas as pd

from ..base.log import print_file_info
from ..base.utils import Timer
from ..calculation.calculate_para import cal_Ep
from ..calculation.calculate_partition_func import cal_Q_nlte_2T


def _exomolhr_meta_filepath():
    """Return the fixed ExoMolHR metadata filepath."""
    return os.path.join(
        os.path.dirname(__file__),
        'meta',
        'ExoMolHR_list.csv',
    )


def _split_meta_list(value):
    """Split a comma-separated metadata field into a clean list."""
    if pd.isna(value):
        return []
    return [item.strip() for item in str(value).split(',') if item.strip() != '']


def _exomolhr_meta_row(molecule, isotopologue):
    """Read the ExoMolHR metadata row for a molecule/isotopologue pair."""
    meta_df = pd.read_csv(_exomolhr_meta_filepath())
    meta_row = meta_df[
        (meta_df['molecule'] == molecule) & (meta_df['iso-slug'] == isotopologue)
    ]
    if meta_row.empty:
        # Fallback to just iso-slug if molecule match fails
        meta_row = meta_df[meta_df['iso-slug'] == isotopologue]
        
    if meta_row.empty:
        raise ValueError(
            f"No ExoMolHR metadata found for molecule={molecule}, isotopologue={isotopologue}."
        )
    return meta_row.iloc[0]


def exomolhr_metadata(molecule, isotopologue):
    """
    Return parsed ExoMolHR metadata for a molecule/isotopologue pair.

    Returns
    -------
    dict
        Metadata including dataset, mass, abundance, column names and formats.
    """
    meta_row = _exomolhr_meta_row(molecule, isotopologue)
    main_col = _split_meta_list(meta_row['Main column'])
    main_fmt = _split_meta_list(meta_row['Main format'])
    qn_label = _split_meta_list(meta_row['QN label'])
    qn_format = _split_meta_list(meta_row['QN format'])

    exomolhr_col = main_col + [label + "'" for label in qn_label] + [label + '"' for label in qn_label]
    exomolhr_fmt = main_fmt + qn_format + qn_format

    return {
        'dataset': str(meta_row['dataset']),
        'mass': float(meta_row['mass']),
        'abundance': float(meta_row['abundance']),
        'qnslabel_list': qn_label,
        'qnsformat_list': qn_format,
        'exomolhr_col': exomolhr_col,
        'exomolhr_fmt': exomolhr_fmt,
    }


def _match_single_dataset_file(filepaths, prefix, suffix, label):
    """Return the single dataset filepath and dataset label."""
    datasets = []
    for filepath in filepaths:
        filename = os.path.basename(filepath)
        if prefix in filename and filename.endswith(suffix):
            # Safer dataset extraction by splitting by '__'
            parts = filename.replace(suffix, '').split('__')
            dataset_name = parts[-1] if len(parts) >= 1 else 'default'
            datasets.append((filepath, dataset_name))
    if len(datasets) == 1:
        return datasets[0]
    if len(datasets) == 0:
        raise FileNotFoundError(f"No ExoMolHR {label} file found.")
    raise ValueError(f"Multiple ExoMolHR {label} files found: {[item[0] for item in datasets]}")


def resolve_exomolhr_filepaths(read_path, molecule, isotopologue):
    """
    Resolve ExoMolHR CSV and PF filepaths using only molecule and isotopologue.
    """
    iso_path = os.path.join(read_path, molecule, isotopologue)
    csv_match = _match_single_dataset_file(
        glob.glob(os.path.join(iso_path, '*.csv')),
        isotopologue,
        '.csv',
        'CSV',
    )
    pf_match = _match_single_dataset_file(
        glob.glob(os.path.join(iso_path, '*.pf')),
        isotopologue,
        '.pf',
        'partition function',
    )
    if csv_match[1] != pf_match[1]:
        # If they don't match, we still allow it but prefer the PF dataset name
        # as it usually follows the official XABC/YNAT naming.
        pass
    return csv_match[0], pf_match[0], pf_match[1]


def resolve_exomolhr_data_info(read_path, molecule, isotopologue):
    """Resolve ``[molecule, isotopologue, dataset]`` for ExoMolHR."""
    _, _, dataset = resolve_exomolhr_filepaths(read_path, molecule, isotopologue)
    return [molecule, isotopologue, dataset]


def read_exomolhr_df(read_path, data_info, min_wn, max_wn, unc_filter):
    """
    Read an ExoMolHR CSV line list into a DataFrame.

    The returned line list uses ``v`` as the transition wavenumber column and
    adds ``E'`` derived from ``E"`` and ``v``.
    """
    t = Timer()
    t.start()
    print('Reading line list ...')

    csv_path, _, _ = resolve_exomolhr_filepaths(read_path, data_info[0], data_info[1])
    meta = exomolhr_metadata(data_info[0], data_info[1])
    exomolhr_df = pd.read_csv(csv_path, header=0, low_memory=False)
    exomolhr_df = exomolhr_df.rename(columns={'nu': 'v'})

    numeric_cols = ['v', 'unc', 'A', 'S', 'E"', "g'", 'g"', "J'", 'J"']
    for col in numeric_cols:
        if col in exomolhr_df.columns:
            exomolhr_df[col] = pd.to_numeric(exomolhr_df[col], errors='coerce')

    exomolhr_df = exomolhr_df[exomolhr_df['v'].between(min_wn, max_wn)]
    if unc_filter != 'None':
        if 'unc' not in exomolhr_df.columns:
            raise ValueError("No uncertainties in ExoMolHR line list. Please do not use uncertainty filter.")
        exomolhr_df = exomolhr_df[exomolhr_df['unc'] <= unc_filter]
    exomolhr_df["E'"] = cal_Ep(exomolhr_df['E"'].values, exomolhr_df['v'].values)
    exomolhr_df = exomolhr_df.reset_index(drop=True)

    t.end()
    print_file_info('Line list', meta['exomolhr_col'], meta['exomolhr_fmt'])
    print('Finished reading line list!\n')
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')
    return exomolhr_df


def read_exomolhr_pf(read_path, data_info, T_list):
    """Read ExoMolHR partition-function values for selected temperatures."""
    _, pf_path, _ = resolve_exomolhr_filepaths(read_path, data_info[0], data_info[1])
    pf_df = pd.read_csv(pf_path, sep=r'\s+', names=['T', 'Q'], header=None, engine='python')
    Q_list = pf_df[pf_df['T'].isin(T_list)]['Q']
    if Q_list.empty:
        raise ValueError(
            'No specified temperature dependent partition funtion value.',
            'Please change the temperature(s) or calculate the partition function at first.',
        )
    return Q_list.to_numpy(dtype=float)


def process_exomolhr_linelist(exomolhr_df):
    """Extract line-list arrays needed for stick-spectra and cross-section calculations."""
    A = exomolhr_df['A'].values
    v = exomolhr_df['v'].values
    Ep = exomolhr_df["E'"].values
    Epp = exomolhr_df['E"'].values
    gp = exomolhr_df["g'"].values
    return A, v, Ep, Epp, gp


def process_exomolhr_linelist_Q(exomolhr_df, T_list, Tvib_list, Trot_list):
    """
    Prepare ExoMolHR line-list arrays and partition functions.

    Supports LTE and two-temperature non-LTE calculations.
    """
    from pyexocross.core import NLTEMethod, vib_label, rot_label, abs_emi, read_path, data_info

    A, v, Ep, Epp, gp = process_exomolhr_linelist(exomolhr_df)

    if NLTEMethod == 'L':
        Q_arr = read_exomolhr_pf(read_path, data_info, T_list)
        Evibp = None
        Erotp = None
        Evibpp = None
        Erotpp = None
    elif NLTEMethod == 'T':
        if len(vib_label) == 0 or len(rot_label) == 0:
            raise ValueError("vib_label and rot_label must be specified for non-LTE calculations with NLTEMethod='T'.")

        all_states_list = []

        vib_cols_p = [label + "'" for label in vib_label]
        rot_cols_p = [label + "'" for label in rot_label]
        upper_states_df = exomolhr_df[vib_cols_p + rot_cols_p + ["E'", "g'"]].drop_duplicates()
        min_rot_energy_upper = upper_states_df.groupby(vib_cols_p)["E'"].min()

        def get_evib_upper(row):
            key = tuple(row[vib_cols_p].values) if len(vib_cols_p) > 1 else row[vib_cols_p[0]]
            return min_rot_energy_upper[key]

        upper_states_df['Evib'] = upper_states_df.apply(get_evib_upper, axis=1)
        upper_states_df['Erot'] = upper_states_df["E'"] - upper_states_df['Evib']
        all_states_list.append(upper_states_df[['Evib', 'Erot', "g'"]].rename(columns={"g'": 'g'}))
        if abs_emi == 'Em':
            Evibp = exomolhr_df.apply(get_evib_upper, axis=1).values
            Erotp = Ep - Evibp
        else:
            Evibp = None
            Erotp = None

        vib_cols_pp = [label + '"' for label in vib_label]
        rot_cols_pp = [label + '"' for label in rot_label]
        lower_states_df = exomolhr_df[vib_cols_pp + rot_cols_pp + ['E"', 'g"']].drop_duplicates()
        min_rot_energy_lower = lower_states_df.groupby(vib_cols_pp)['E"'].min()

        def get_evib_lower(row):
            key = tuple(row[vib_cols_pp].values) if len(vib_cols_pp) > 1 else row[vib_cols_pp[0]]
            return min_rot_energy_lower[key]

        lower_states_df['Evib'] = lower_states_df.apply(get_evib_lower, axis=1)
        lower_states_df['Erot'] = lower_states_df['E"'] - lower_states_df['Evib']
        all_states_list.append(lower_states_df[['Evib', 'Erot', 'g"']].rename(columns={'g"': 'g'}))
        if abs_emi == 'Ab':
            Evibpp = exomolhr_df.apply(get_evib_lower, axis=1).values
            Erotpp = Epp - Evibpp
        else:
            Evibpp = None
            Erotpp = None

        all_states_df = pd.concat(all_states_list, ignore_index=True)
        Q_arr = cal_Q_nlte_2T(Tvib_list, Trot_list, all_states_df['Evib'].values, all_states_df['Erot'].values, all_states_df['g'].values)
    else:
        raise ValueError("ExoMolHR line-list calculations support LTE and two-temperature non-LTE (NLTEMethod='T') only.")

    return A, v, Ep, Epp, gp, Q_arr, Evibp, Erotp, Evibpp, Erotpp
