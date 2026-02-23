"""
Load HITRAN database files.

This module provides functions for reading HITRAN parameter files and
processing HITRAN linelist data.
"""
import os
import numpy as np
import pandas as pd
from ..base.utils import Timer
from ..base.log import print_file_info
from ..base.constants import c2, Tref, Pref
from ..calculation.calculate_para import cal_Ep_hitran, cal_Jp, cal_F
from ..calculation.calculate_partition_func import cal_Q_nlte_2T
from ..process.hitran_qn import separate_QN_hitran


def _expected_hitran_basename():
    """Return expected HITRAN base filename: molecule__isoslug."""
    from pyexocross.core import data_info
    molecule, isoslug = data_info[0], data_info[1]
    return f"{molecule}__{isoslug}"


def _resolve_hitran_par_path(read_path):
    """
    Resolve and validate HITRAN .par path.

    Supported input:
    - directory path: must contain <molecule>__<isoslug>.par
    - file path: must be exactly <molecule>__<isoslug>.par
    """
    expected_base = _expected_hitran_basename()
    expected_name = f"{expected_base}.par"
    clean_path = read_path.rstrip('/')

    if os.path.isdir(clean_path):
        # <read_path>/<molecule>/<isoslug>/<molecule>__<isoslug>.par
        from pyexocross.core import data_info
        molecule, isoslug = data_info[0], data_info[1]
        candidates = [
            os.path.join(clean_path, molecule, isoslug, expected_name),
        ]
        for par_path in candidates:
            if os.path.exists(par_path):
                return par_path
        raise ValueError(
            "Missing HITRAN line list file. Expected one of: "
            f"{candidates[0]} or {candidates[1]}"
        )

    if not os.path.exists(clean_path):
        raise ValueError('The input file ' + clean_path + ' does not exist.')

    fname = os.path.basename(clean_path)
    if fname != expected_name:
        raise ValueError(
            f"Invalid HITRAN line list filename: {fname}. "
            f"Expected: {expected_name}"
        )
    return clean_path


def _resolve_hitran_pf_path(read_path):
    """
    Resolve and validate HITRAN partition-function path.

    Expected filename (same base as line list):
    - <molecule>__<isoslug>.pf
    - <molecule>__<isoslug>.txt
    """
    par_path = _resolve_hitran_par_path(read_path)
    base, _ = os.path.splitext(par_path)
    candidates = [base + '.pf', base + '.txt']
    for pf_path in candidates:
        if os.path.exists(pf_path):
            return pf_path
    raise ValueError(
        "No partition function file found. Expected one of: "
        f"{candidates[0]} or {candidates[1]}"
    )


### Read HITRAN Linelist File
def read_parfile(read_path):
    """
    Read HITRAN .par file in chunks.

    Reads the HITRAN parameter file in chunks for memory efficiency.

    Parameters
    ----------
    read_path : str
        HITRAN read path. Can be either:
        - folder containing <molecule>__<isoslug>.par
        - exact path to <molecule>__<isoslug>.par

    Returns
    -------
    pd.DataFrame
        Complete HITRAN parameter file DataFrame
    """
    par_path = _resolve_hitran_par_path(read_path)
    # Initialise the iterator object.
    read_par = pd.read_csv(par_path, chunksize=100_000, iterator=True, header=None, encoding='utf-8')
    par_df = pd.DataFrame()
    for chunk in read_par:
        par_df = pd.concat([par_df, chunk])
    return par_df

# Process HITRAN linelist data
def read_hitran_parfile(read_path, parfile_df, minv, maxv, unclimit, Slimit):
    """
    Parse HITRAN2004 format parameter file and extract molecular line data.

    Reads and parses fixed-format HITRAN2004 .par file, extracts all columns,
    and filters by molecule ID, isotopologue ID, wavenumber range, uncertainty,
    and intensity threshold.

    Parameters
    ----------
    read_path : str
        HITRAN read path (folder or exact .par path)
    parfile_df : pd.DataFrame
        Raw parameter file DataFrame (single column with fixed-width strings)
    minv : float
        Minimum wavenumber in cm^-1
    maxv : float
        Maximum wavenumber in cm^-1
    unclimit : float or str
        Uncertainty limit ('None' to disable)
    Slimit : float or str
        Intensity threshold ('None' to disable)

    Returns
    -------
    pd.DataFrame
        Parsed HITRAN DataFrame with columns: M, I, v, S, A, gamma_air, gamma_self,
        Epp, n_air, delta_air, Vp, Vpp, Qp, Qpp, Unc, gp, gpp
    """    
    from pyexocross.core import species_main_id, species_sub_id
    print('Reading HITRAN2004 format data file ...')
    t = Timer()
    t.start()
    par_path = _resolve_hitran_par_path(read_path)
    par_filename = os.path.basename(par_path)
    if (len(str(parfile_df[0][0])) < 160):
        raise ValueError('The file ' + par_filename + ' is not a HITRAN2004 format data file.')
    hitran_column_name = ['M','I','v','S','A','gamma_air','gamma_self',
                          'Epp','n_air','delta_air','Vp','Vpp','Qp','Qpp',
                          'Ierr','Iref','flag','gp','gpp']

    hitran_df = pd.DataFrame()
    hitran_df['M'] = pd.to_numeric(parfile_df[0].map(lambda x: x[0:2]), errors='coerce').astype('int64')                 # Molecule identification number
    hitran_df['I'] = pd.to_numeric(parfile_df[0].map(lambda x: x[2:3]), errors='coerce').astype('int64')                 # Isotopologue number
    hitran_df['v'] = pd.to_numeric(parfile_df[0].map(lambda x: x[3:15]), errors='coerce').astype('float64')              # Transition wavenumber (in cm^{-1})
    hitran_df['S'] = pd.to_numeric(parfile_df[0].map(lambda x: x[15:25]), errors='coerce').astype('float64')             # Intensity (cm^{-1} / (molecule cm^{-2}))
    hitran_df['A'] = pd.to_numeric(parfile_df[0].map(lambda x: x[25:35]), errors='coerce').astype('float64')             # The Einstein-A coefficient (s^{-1}) of a transition
    hitran_df['gamma_air'] = pd.to_numeric(parfile_df[0].map(lambda x: x[35:40]), errors='coerce').astype('float64')     # Air-broadened half-width at half maximum (HWHM) coefficient (cm^{-1} atm^{-1})
    hitran_df['gamma_self'] = pd.to_numeric(parfile_df[0].map(lambda x: x[40:45]), errors='coerce').astype('float64')    # Self-broadened half-width at half maximum (HWHM) coefficient (cm^{-1} atm^{-1})
    hitran_df['Epp'] = pd.to_numeric(parfile_df[0].map(lambda x: x[45:55]), errors='coerce').astype('float64')           # Lower state energy (cm^{-1})
    hitran_df['n_air'] = pd.to_numeric(parfile_df[0].map(lambda x: x[55:59]), errors='coerce').astype('float64')         # Temperature-dependent exponent for gamma_air
    hitran_df['delta_air'] = pd.to_numeric(parfile_df[0].map(lambda x: x[59:67]), errors='coerce').astype('float64')     # Air pressure_include line shift (cm^{-1} atm^{-1})
    hitran_df['Vp'] = parfile_df[0].map(lambda x: x[67:82])                                                              # Upper-state "global" quanta
    hitran_df['Vpp'] = parfile_df[0].map(lambda x: x[82:97])                                                             # Lower-state "global" quanta
    hitran_df['Qp'] = parfile_df[0].map(lambda x: x[97:112])                                                             # Upper-state "local" quanta
    hitran_df['Qpp'] = parfile_df[0].map(lambda x: x[112:127])                                                           # Lower-state "local" quanta
    #hitran_df['Unc'] = parfile_df[0].map(lambda x: x[127:128])                                                          # Uncertainty code, first integer in the error code
    hitran_df['Unc'] = pd.to_numeric(parfile_df[0].map(lambda x: x[127:128]), errors='coerce').astype('int64')           # Uncertainty code, first integer in the error code
    hitran_df['gp'] = pd.to_numeric(parfile_df[0].map(lambda x: x[146:153]), errors='coerce').astype('int64')            # Statistical weight of the upper state
    hitran_df['gpp'] = pd.to_numeric(parfile_df[0].map(lambda x: x[153:160]), errors='coerce').astype('int64')           # Statistical weight of the lower state
    
    hitran_df = hitran_df[hitran_df['M'].isin([species_main_id])]
    hitran_df = hitran_df[hitran_df['I'].isin([species_sub_id])]
    hitran_df = hitran_df[hitran_df['v'].between(minv, maxv)]
    if unclimit != 'None':
        hitran_df = hitran_df[hitran_df['Unc'] >= int(('%e' % unclimit)[-1])]
    if Slimit != 'None':
        hitran_df = hitran_df[hitran_df['S'] >= Slimit]
    t.end()
    hitran_col = ['M','I','v','S','A','gamma_air','gamma_self','E"','n_air','delta_air',
                  "V'",'V"',"Q'",'Q"','Ierr','Iref','flag',"g'",'g"']
    hitran_fmt = ['%2d','%1d','%12.6f','%10.3e','%10.3e','%5.4f','%5.4f','%10.4f','%4.2f','%8.6f',
                  '%15s','%15s','%15s','%15s','%6d','%12d','%1s','%7.1f','%7.1f']
    print_file_info('HITRAN2004 format par', hitran_col, hitran_fmt)
    return hitran_df

# ### Read Partition Function File From HITRANOnline
# # Read partition function file from HITRANOnline
# def read_hitranweb_pf(T_list):
#     isometa_url = 'https://hitran.org/docs/iso-meta/'
#     iso_meta_table = pd.read_html(isometa_url)[species_main_id - 1]
#     iso_meta_row = iso_meta_table[iso_meta_table['local ID'].isin([species_sub_id])]
#     #Q_ref = float(iso_meta_row.loc[0][6].replace('\xa0×\xa010','E'))
#     Q_url = 'https://hitran.org/data/Q/' + iso_meta_row.loc[0][7]
#     Q_content = requests.get(Q_url).text
#     Q_col_name = ['T', 'Q']
#     pf_df = pd.read_csv(StringIO(Q_content), sep='\\s+', names=Q_col_name, header=None)
#     try:
#         Q_list = pf_df[pf_df['T'].isin(T_list)]['Q']
#     except:
#         raise ValueError('No specified temperature dependent partition funtion value.', 
#                           'Please change the temperature(s) or calculate the partition function at first.')
#     Q_arr = Q_list.to_numpy(dtype=float)  
#     return Q_arr

# Read partition function file from HITRANOnline local file
def read_hitran_pf(T_list):
    """
    Read HITRAN partition function file and return Q(T) for specified temperatures.

    Parameters
    ----------
    T_list : list of float
        List of temperatures in Kelvin

    Returns
    -------
    np.ndarray
        Partition function values Q(T) for each temperature, shape (n_temps,)
    """
    from pyexocross.core import read_path
    Q_col_name = ['T', 'Q']
    pf_path = _resolve_hitran_pf_path(read_path)
    pf_df = pd.read_csv(pf_path, sep='\\s+', names=Q_col_name, header=None)
    try:
        Q_list = pf_df[pf_df['T'].isin(T_list)]['Q']
    except:
        raise ValueError('No specified temperature dependent partition funtion value.', 
                          'Please change the temperature(s) or calculate the partition function at first.')
    Q_arr = Q_list.to_numpy(dtype=float)    
    return Q_arr
    
def process_hitran_linelist(hitran_linelist_df):
    """
    Extract line parameters from HITRAN linelist DataFrame.

    Extracts transition parameters needed for line profile calculations.

    Parameters
    ----------
    hitran_linelist_df : pd.DataFrame
        HITRAN linelist DataFrame with processed columns

    Returns
    -------
    tuple
        A tuple containing:
        - A : np.ndarray, Einstein A coefficients
        - v : np.ndarray, transition wavenumbers
        - Ep : np.ndarray, upper state energies
        - Epp : np.ndarray, lower state energies
        - gp : np.ndarray, upper state degeneracies
        - n_air : np.ndarray, temperature exponents
        - gamma_air : np.ndarray, air-broadening coefficients
        - gamma_self : np.ndarray, self-broadening coefficients
        - delta_air : np.ndarray, air pressure shift coefficients
    """
    A = hitran_linelist_df['A'].values
    Ep = hitran_linelist_df["E'"].values
    Epp = hitran_linelist_df['E"'].values
    n_air = hitran_linelist_df['n_air'].values
    gamma_air = hitran_linelist_df['gamma_air'].values
    gamma_self = hitran_linelist_df['gamma_self'].values
    delta_air = hitran_linelist_df['delta_air'].values
    gp = hitran_linelist_df["g'"].values
    v = hitran_linelist_df['v'].values
    # if broad == 'Air':
    #   v = hitran_df['v'].values + delta_air * (P - P_ref) / P
    # else:
    #   v = hitran_df['v'].values
    return (A, v, Ep, Epp, gp, n_air, gamma_air, gamma_self, delta_air)

# Prepare HITRAN data for calculations (extract Evib/Erot, calculate partition functions)
def process_hitran_linelist_Q(hitran_linelist_df, T_list, Tvib_list, Trot_list):
    """
    Prepare HITRAN data for calculations: extract Evib/Erot and calculate partition functions.

    This function is called once before parallel calculations to avoid repeating
    data processing. Handles both LTE and non-LTE (two-temperature) methods.

    Parameters
    ----------
    hitran_linelist_df : pd.DataFrame
        HITRAN linelist DataFrame with processed columns
    T_list : list of float
        List of temperatures for LTE calculations (empty for non-LTE)
    Tvib_list : list of float
        List of vibrational temperatures for non-LTE (empty [] for LTE)
    Trot_list : list of float
        List of rotational temperatures for non-LTE (empty [] for LTE)

    Returns
    -------
    tuple
        A tuple containing:
        - A : np.ndarray, Einstein A coefficients
        - v : np.ndarray, transition wavenumbers
        - Ep : np.ndarray, upper state energies
        - Epp : np.ndarray, lower state energies
        - gp : np.ndarray, upper state degeneracies
        - n_air : np.ndarray, temperature exponents
        - gamma_air : np.ndarray, air-broadening coefficients
        - gamma_self : np.ndarray, self-broadening coefficients
        - delta_air : np.ndarray, air pressure shift coefficients
        - Q_arr : np.ndarray, partition function array
        - Evibp : np.ndarray or None, upper state vibrational energies
        - Erotp : np.ndarray or None, upper state rotational energies
        - Evibpp : np.ndarray or None, lower state vibrational energies
        - Erotpp : np.ndarray or None, lower state rotational energies
    """
    from pyexocross.core import NLTEMethod, vib_label, rot_label, abs_emi
    # Extract line data ONCE
    A, v, Ep, Epp, gp, n_air, gamma_air, gamma_self, delta_air = process_hitran_linelist(hitran_linelist_df)

    # Calculate partition functions
    if NLTEMethod == 'L':
        # LTE: Tvib_list and Trot_list are empty lists []
        Q_arr = read_hitran_pf(T_list)
        Evibp = None
        Erotp = None
        Evibpp = None
        Erotpp = None
    elif NLTEMethod == 'T':
        # Non-LTE T mode: only use Tvib_list and Trot_list, not T_list
        if len(vib_label) == 0 or len(rot_label) == 0:
            raise ValueError("vib_label and rot_label must be specified for non-LTE calculations with NLTEMethod='T'")
        
        # Extract Evib and Erot using ExoMol approach (similar to lines 2536-2547)
        # For HITRAN, we need to handle both upper and lower states separately
        
        # Collect all unique states (both upper and lower) for Qnlte calculation
        all_states_list = []
        
        # Process upper states
        if all(col in hitran_linelist_df.columns for col in [v + "'" for v in vib_label] + [r + "'" for r in rot_label] + ["E'", "g'"]):
            vib_cols_p = [v + "'" for v in vib_label]
            rot_cols_p = [r + "'" for r in rot_label]
            # Get unique upper states
            upper_states_df = hitran_linelist_df[vib_cols_p + rot_cols_p + ["E'", "g'"]].drop_duplicates()
            # Find minimum energy for each vibrational level (ExoMol approach: groupby vib_label, find min E)
            min_rot_energy_upper = upper_states_df.groupby(vib_cols_p)["E'"].min()
            # Calculate Evib and Erot for upper states
            def get_evib_upper(row):
                key = tuple(row[vib_cols_p].values) if len(vib_cols_p) > 1 else row[vib_cols_p[0]]
                return min_rot_energy_upper[key]
            upper_states_df['Evib'] = upper_states_df.apply(get_evib_upper, axis=1)
            upper_states_df['Erot'] = upper_states_df["E'"] - upper_states_df['Evib']
            all_states_list.append(upper_states_df[['Evib', 'Erot', "g'"]].rename(columns={"g'": 'g'}))
            
            # Map Evib and Erot back to transitions for upper states
            if abs_emi == 'Em':
                def get_evib_trans_upper(row):
                    key = tuple(row[vib_cols_p].values) if len(vib_cols_p) > 1 else row[vib_cols_p[0]]
                    return min_rot_energy_upper[key]
                Evibp = hitran_linelist_df.apply(get_evib_trans_upper, axis=1).values
                Erotp = Ep - Evibp
            else:
                Evibp = None
                Erotp = None
        else:
            Evibp = None
            Erotp = None
        
        # Process lower states
        if all(col in hitran_linelist_df.columns for col in [v + '"' for v in vib_label] + [r + '"' for r in rot_label] + ['E"', 'g"']):
            vib_cols_pp = [v + '"' for v in vib_label]
            rot_cols_pp = [r + '"' for r in rot_label]
            # Get unique lower states
            lower_states_df = hitran_linelist_df[vib_cols_pp + rot_cols_pp + ['E"', 'g"']].drop_duplicates()
            # Find minimum energy for each vibrational level (ExoMol approach)
            min_rot_energy_lower = lower_states_df.groupby(vib_cols_pp)['E"'].min()
            # Calculate Evib and Erot for lower states
            def get_evib_lower(row):
                key = tuple(row[vib_cols_pp].values) if len(vib_cols_pp) > 1 else row[vib_cols_pp[0]]
                return min_rot_energy_lower[key]
            lower_states_df['Evib'] = lower_states_df.apply(get_evib_lower, axis=1)
            lower_states_df['Erot'] = lower_states_df['E"'] - lower_states_df['Evib']
            all_states_list.append(lower_states_df[['Evib', 'Erot', 'g"']].rename(columns={'g"': 'g'}))
            
            # Map Evib and Erot back to transitions for lower states
            if abs_emi == 'Ab':
                def get_evib_trans_lower(row):
                    key = tuple(row[vib_cols_pp].values) if len(vib_cols_pp) > 1 else row[vib_cols_pp[0]]
                    return min_rot_energy_lower[key]
                Evibpp = hitran_linelist_df.apply(get_evib_trans_lower, axis=1).values
                Erotpp = Epp - Evibpp
            else:
                Evibpp = None
                Erotpp = None
        else:
            Evibpp = None
            Erotpp = None
        
        # Combine all unique states for Qnlte calculation
        if len(all_states_list) > 0:
            all_states_df = pd.concat(all_states_list, ignore_index=True)
            all_states_Evib = all_states_df['Evib'].values
            all_states_Erot = all_states_df['Erot'].values
            all_states_g = all_states_df['g'].values
            # Calculate NLTE partition function (unified Q_arr)
            Q_arr = cal_Q_nlte_2T(Tvib_list, Trot_list, all_states_Evib, all_states_Erot, all_states_g)
        else:
            raise ValueError("Could not extract unique states for NLTE partition function calculation")
    else:
        raise ValueError("Please choose: 'LTE' or 'Non-LTE' and 'T' two temperatures (Tvib, Trot) method.")
    
    return (A, v, Ep, Epp, gp, n_air, gamma_air, gamma_self, delta_air, 
            Q_arr, Evibp, Erotp, Evibpp, Erotpp)
