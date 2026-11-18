"""
Load ExoMol database files.

This module provides functions for reading ExoMol states files, transition files,
partition function files, and broadening parameter files.
"""
import re
import os
import glob
import numpy as np
import pandas as pd
import dask.dataframe as dd
from ..base.utils import Timer
from ..base.log import print_file_info
from ..base.constants import LARGE_TRANS_FILE_BYTES
from ..base.large_file import command_decompress
from ..base.qn_metadata import normalized_states_columns
from ..process.exomol_exomolhr_broad import (
    broad_required_line_columns,
    read_broad_file,
    read_broad,
    extract_broad,
)

# Read Input Files
## Read ExoMol Database Files
### Read States File
def preferred_files(folder, suffix):
    """Return matching files with uncompressed variants preferred over .bz2."""
    selected = {}
    plain_files = glob.glob(os.path.join(folder, f'*{suffix}'))
    compressed_files = glob.glob(os.path.join(folder, f'*{suffix}.bz2'))
    for filepath in plain_files + compressed_files:
        logical_path = filepath[:-4] if filepath.endswith('.bz2') else filepath
        if logical_path not in selected or not filepath.endswith('.bz2'):
            selected[logical_path] = filepath
    return sorted(selected.values())


def get_statesfile(read_path, data_info):
    """Return the preferred states file, choosing .states over .states.bz2."""
    folder = read_path + '/'.join(data_info) + '/'
    basename = '__'.join(data_info[-2:]) + '.states'
    plain = os.path.join(folder, basename)
    compressed = plain + '.bz2'
    if os.path.exists(plain):
        return plain
    if os.path.exists(compressed):
        return compressed
    raise ValueError("No such states file, please check the read path and states filename format!")


def read_all_states(read_path, data_info, check_uncertainty, states_col, states_fmt):
    """
    Read complete ExoMol states file.

    Reads the states file (compressed or uncompressed) and returns a DataFrame
    with state information including id, energy, degeneracy, J, and optionally
    uncertainty, lifetime, and g-factor.

    Parameters
    ----------
    read_path : str
        Base path to ExoMol database directory

    Returns
    -------
    pd.DataFrame
        States DataFrame with columns: id, E, g, J, and optionally unc, tau, gfac
    """
    t = Timer()
    t.start()
    print('Reading states ...')
    states_filename = get_statesfile(read_path, data_info)
    read_kwargs = dict(
        sep=r'\s+',
        header=None,
        chunksize=500_000,
        iterator=True,
        low_memory=False,
        dtype=object,
    )
    if states_filename.endswith('.bz2'):
        read_kwargs['compression'] = 'bz2'
    states_df = pd.concat(pd.read_csv(states_filename, **read_kwargs), ignore_index=True)
    if check_uncertainty == True:
        states_df = states_df.rename(columns={0:'id',1:'E',2:'g',3:'J',4:'unc'})
        # Use float64 for Unc column to avoid precision issues with very small values
        convert_dict = {'id':np.int32,'E':np.float64,'g':np.int32,'J':np.float16,'unc':np.float64}
    else:      
        states_df = states_df.rename(columns={0:'id',1:'E',2:'g',3:'J'})  
        convert_dict = {'id':np.int32,'E':np.float64,'g':np.int32,'J':np.float16}
    states_df = states_df[states_df['E'].notna()]
    states_df = states_df.astype(convert_dict)
    t.end('load')
    print_file_info('States', states_col, states_fmt)
    print('Finished reading states!\n')       
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')                
    return states_df

# Read part states file for calculating LTE or non-LTE
def read_part_states(
    states_df,
    UncFilter,
    NLTEMethod,
    NLTEPath,
    check_uncertainty,
    check_lifetime,
    check_gfactor,
    QNslabel_list,
    states_col,
    states_fmt,
    QNs_label,
    QNsFilter,
    QNs_value,
    vib_label,
    rot_label,
):
    """
    Read and filter states DataFrame for LTE or non-LTE calculations.

    Applies uncertainty filtering, selects relevant columns, and processes
    states according to the specified non-LTE method (LTE, two-temperature,
    custom density, or custom population).

    Parameters
    ----------
    states_df : pd.DataFrame
        Complete states DataFrame from read_all_states()

    Returns
    -------
    pd.DataFrame
        Filtered and processed states DataFrame ready for calculations
    """
    colname = normalized_states_columns(states_col, states_df.shape[1]) if states_col else list(states_df.columns)
    available_cols = set(colname)
    col_unc = ['unc'] if check_uncertainty == True and 'unc' in available_cols else []
    col_lifetime = ['tau'] if check_lifetime == True and 'tau' in available_cols else []
    col_gfac = ['gfac'] if check_gfactor == True and 'gfac' in available_cols else []

    states_parts_df = states_df.loc[:, states_df.columns[:len(colname)]].copy()
    states_parts_df.columns = colname
    if UncFilter != 'None':
        if col_unc:
            states_parts_df = states_parts_df[states_parts_df['unc'] <= UncFilter]
            # states_parts_df.set_index(['id'], inplace=True, drop=False)
            # states_parts_df['id'] = pd.to_numeric(states_parts_df['id'])
        else:
            raise ValueError("No uncertainties in states file. Please do not use uncertainty filter.")

    if not states_col:
        colname = ['id','E','g','J'] + col_unc + col_lifetime + col_gfac + QNslabel_list
        states_parts_df = states_parts_df.loc[:, states_parts_df.columns[:len(colname)]].copy()
        states_parts_df.columns = colname
        available_cols = set(colname)
    optional_cols = col_unc + col_lifetime + col_gfac
    # Keep only the required columns and avoid SettingWithCopyWarning by working on a copy
    QNcolumns = [
        col
        for col in list(dict.fromkeys(
            ['id', 'E', 'g', 'J'] + optional_cols + QNslabel_list + QNs_label
        ))
        if col in states_parts_df.columns
    ]

    # LTE 
    if NLTEMethod == 'L':
        states_part_df = states_parts_df[QNcolumns]
    # Non-LTE using two temperatures
    elif NLTEMethod == 'T':
        QNcolumns_2Ts = ['id','E','g','J'] + optional_cols + QNs_label + rot_label + vib_label
        # Remove duplicates while preserving order
        QNcolumns_2T = [col for col in list(dict.fromkeys(QNcolumns_2Ts)) if col in available_cols]
        states_part_df = states_parts_df[QNcolumns_2T].copy()
        # Define the vibronic origin independently for every vibrational/
        # electronic manifold. Different TiO electronic states can have
        # different minimum J and parity, so using the global ground-state
        # rotational labels incorrectly removes most manifolds.
        states_part_df['Evib'] = states_part_df.groupby(
            vib_label,
            dropna=False,
        )['E'].transform('min')
        states_part_df['Erot'] = states_part_df['E'] - states_part_df['Evib']
        output_columns = [col for col in QNcolumns if col in states_part_df.columns]
        states_part_df = states_part_df[output_columns + ['Erot', 'Evib']]
    # Non-LTE using custom density
    elif NLTEMethod == 'D':
        # Read density file
        nlte_density = pd.read_csv(NLTEPath, header=None, names=['id','density'])
        nlte_density['nvib'] = nlte_density['density'] / nlte_density['density'].sum()
        nlte_density = nlte_density[['id','nvib']]
        states_part_df = pd.merge(states_parts_df, nlte_density, on='id', how='inner')
    # Non-LTE using custom population   
    elif NLTEMethod == 'P':
        # Read population file
        nlte_population = pd.read_csv(NLTEPath, header=None, names=['id','population'])
        nlte_population['pop'] = nlte_population['population'] / nlte_population['population'].sum()
        nlte_population = nlte_population[['id','pop']]
        states_part_df = pd.merge(states_parts_df, nlte_population, on='id', how='inner')
    else:
        raise ValueError("Please choose one non-LTE method from: 'T', 'D' or 'P'.")   
    states_parts_df.set_index(['id'], inplace=True, drop=False) 
    # Use quantum number filters
    if QNsFilter !=[]:    
        for i in range(len(QNs_label)):
            if QNs_value[i] != ['']:
                list_QN_value = str(QNs_value[i]).replace("', '",",").replace("'","").replace('[','').replace(']','').split(',')
                if '' in list_QN_value:
                    states_part_df = states_part_df
                else:
                    states_part_df = states_part_df[states_part_df[QNs_label[i]].isin(list_QN_value)]
    else:
        pass
    states_part_df.index.name = 'index'
    # Work on a copy to avoid SettingWithCopyWarning when assigning to columns
    states_part_df = states_part_df.copy()
    try:
        states_part_df['tau'] = states_part_df['tau'].astype('float')
        states_part_df['gfac'] = states_part_df['gfac'].astype('float')
    except Exception:
        pass
    return states_part_df

### Get transitions File
def get_transfiles(read_path, data_info, prepare=True):
    """
    Get all transition file paths from ExoMol database directory.

    Finds all transition files matching the dataset, filters out old version files,
    and handles decompression for large files.

    Parameters
    ----------
    read_path : str
        Base path to ExoMol database directory
    data_info : list[str]
        Dataset info used to build transition file paths

    Returns
    -------
    list of str
        List of transition file paths (decompressed if necessary)
    """
    # Get all the transitions files from the folder including the older version files which are named by vn(version number).
    print()
    print('Search transitions files ...')
    folder = read_path + '/'.join(data_info) + '/'
    trans_filepaths_all = preferred_files(folder, '.trans')
    num_transfiles_all = len(trans_filepaths_all)    # The number of all transitions files including the older version files.
    trans_filepaths = []    # The list of the lastest transitions files.
    all_decompress_num = 0
    decompress_num = 0
    for i in range(num_transfiles_all):
        split_version = trans_filepaths_all[i].split('__')[-1].split('.')[0]        # Split the filenames.
        # There are five format filenames.
        # The lastest transitions files named in five formats:
        # 1. Filenames are named with the name of isotopologue and dataset. 
        #    End with .trans.bz2.
        #    e.g. 14N-16O__XABC.trans.bz2'
        # 2. Filenames are named with the name of isotopologue and dataset. 
        #    Also have the range of wavenumbers xxxxx-yyyyy.
        #    End with .trans.bz2.
        #    e.g. 1H2-16O__POKAZATEL__00000-00100.trans.bz2
        # 3. Filenames are named with the name of isotopologue and dataset. 
        #    Also have E2, M1, etc after dataset name. 
        #    e.g. 16O2__SWYT__E2.trans.bz2
        # 4. The older version transitions files are named with vn(version number) based on the first format of the lastest files.
        #    e.g. 14N-16O__XABC_v2.trans.bz2
        # 5. The older version transitions files are named with updated date (yyyymmdd).
        #    e.g. 1H3_p__MiZATeP__20170330.trans.bz2
        # After split the filenames:
        # The first format filenames only leave the dataset name, e.g. XABC.
        # The second format filenames only leave the range of the wavenumber, e.g. 00000-00100.
        # The third format filenames leave two parts(dataset name and E2, M1, etc), e.g. SWYT and E2.
        # The fourth format filenames leave two parts(dataset name and version number), e.g. XABC and v2.
        # The fifth format filenames only leave the updated date, e.g. 20170330.
        # This program only process the lastest data, so extract the filenames named by the first three formats.
        if (
            split_version == data_info[-1]
            or re.fullmatch(r'\d{5}-\d{5}', split_version)
            or re.fullmatch(r'[A-Z]\d+', split_version)
        ):   
            original_path = trans_filepaths_all[i]
            file_size_bytes = os.path.getsize(original_path)
            trans_filepath = original_path
            if prepare and original_path.endswith('.bz2') and file_size_bytes >= LARGE_TRANS_FILE_BYTES:
                (trans_filepath, num_dec) = command_decompress(original_path)
                all_decompress_num += 1
                decompress_num += num_dec
            trans_filepaths.append(trans_filepath)
        else:
            pass
    print('{:45s} : {}'.format('Number of all transitions files', num_transfiles_all))
    print('{:45s} : {}'.format('Number of all decompressed transitions files', all_decompress_num))
    print('{:45s} : {}'.format('Number of new decompressed transitions files', decompress_num))
    print()
    return trans_filepaths  

def get_part_transfiles(read_path, data_info, min_wn, max_wn, prepare=True):
    """
    Get transition file paths filtered by wavenumber range.

    Finds transition files matching the dataset and wavenumber range,
    filters out old version files, and handles decompression for large files.

    Parameters
    ----------
    read_path : str
        Base path to ExoMol database directory
    data_info : list[str]
        Dataset info used to build transition file paths
    min_wn : float
        Minimum wavenumber in cm^-1 used for selecting transition files
    max_wn : float
        Maximum wavenumber in cm^-1 used for selecting transition files

    Returns
    -------
    list of str
        List of transition file paths within the specified wavenumber range
    """
    # Get all the transitions files from the folder including the older version files which are named by vn(version number).
    print()
    print('Search transitions files in specified range ...')
    folder = read_path + '/'.join(data_info) + '/'
    trans_filepaths_all = preferred_files(folder, '.trans')
    num_transfiles_all = len(trans_filepaths_all)    # The number of all transitions files including the older version files.
    trans_filepaths = []    # The list of the lastest transitions files.
    all_decompress_num = 0
    decompress_num = 0
    for i in range(num_transfiles_all):
        split_version = trans_filepaths_all[i].split('__')[-1].split('.')[0]        # Split the filenames.
        # There are five format filenames.
        # The lastest transitions files named in five formats:
        # 1. Filenames are named with the name of isotopologue and dataset. 
        #    End with .trans.bz2.
        #    e.g. 14N-16O__XABC.trans.bz2'
        # 2. Filenames are named with the name of isotopologue and dataset. 
        #    Also have the range of wavenumbers xxxxx-yyyyy.
        #    End with .trans.bz2.
        #    e.g. 1H2-16O__POKAZATEL__00000-00100.trans.bz2
        # 3. Filenames are named with the name of isotopologue and dataset. 
        #    Also have E2, M1, etc after dataset name. 
        #    e.g. 16O2__SWYT__E2.trans.bz2
        # 4. The older version transitions files are named with vn(version number) based on the first format of the lastest files.
        #    e.g. 14N-16O__XABC_v2.trans.bz2
        # 5. The older version transitions files are named with updated date (yyyymmdd).
        #    e.g. 1H3_p__MiZATeP__20170330.trans.bz2
        # After split the filenames:
        # The first format filenames only leave the dataset name, e.g. XABC.
        # The second format filenames only leave the range of the wavenumber, e.g. 00000-00100.
        # The third format filenames leave two parts(dataset name and E2, M1, etc), e.g. SWYT and E2.
        # The fourth format filenames leave two parts(dataset name and version number), e.g. XABC and v2.
        # The fifth format filenames only leave the updated date, e.g. 20170330.
        # This program only process the lastest data, so extract the filenames named by the first three formats.  
        if split_version == data_info[-1] or re.fullmatch(r'[A-Z]\d+', split_version):
            trans_filepath = trans_filepaths_all[i]
            file_size_bytes = os.path.getsize(trans_filepath)
            if prepare and trans_filepath.endswith('.bz2') and file_size_bytes >= LARGE_TRANS_FILE_BYTES:
                (trans_filepath, num_dec) = command_decompress(trans_filepath)
                all_decompress_num += 1
                decompress_num += num_dec
            trans_filepaths.append(trans_filepath)
        elif re.fullmatch(r'\d{5}-\d{5}', split_version):
            lower = int(split_version.split('-')[0])
            upper = int(split_version.split('-')[1])
            if lower <= max_wn and upper >= min_wn:
                file_size_bytes = os.path.getsize(trans_filepaths_all[i])
                trans_filepath = trans_filepaths_all[i]
                if prepare and trans_filepath.endswith('.bz2') and file_size_bytes >= LARGE_TRANS_FILE_BYTES:
                    (trans_filepath, num_dec) = command_decompress(trans_filepath)
                    all_decompress_num += 1
                    decompress_num += num_dec
                trans_filepaths.append(trans_filepath)
        else:
            pass
    print('{:45s} : {}'.format('Number of all transitions files', num_transfiles_all))
    print('{:45s} : {}'.format('Number of selected transitions files', len(trans_filepaths)))
    print('{:45s} : {}'.format('Number of all decompressed transitions files', all_decompress_num))
    print('{:45s} : {}'.format('Number of new decompressed transitions files', decompress_num))
    print()
    return trans_filepaths

# Read partition function with online webpage
# def read_exomolweb_pf(T_list):
#     pf_url = 'http://www.exomol.com/db/' + '/'.join(data_info) + '/' + '__'.join(data_info[-2:]) + '.pf'
#     pf_col_name = ['T', 'Q']
#     try:
#         pf_content = requests.get(pf_url).text
#         pf_df = pd.read_csv(StringIO(pf_content), sep='\\s+', names=pf_col_name, header=None)
#     except:
#         raise ValueError('No partition function file. Please check the webpage.')
#     try:
#         Q_list = pf_df[pf_df['T'].isin(T_list)]['Q']
#     except:
#         raise ValueError('No specified temperature dependent partition funtion value.', 
#                           'Please change the temperature or calculate the partition function at first.')
#     Q_arr = Q_list.to_numpy(dtype=float)
#     return Q_list 

# Read partition function with local partition function file
def read_exomol_pf(read_path, data_info, T_list, states_df=None):
    """
    Read ExoMol partition function file and return Q(T) for specified temperatures.

    Parameters
    ----------
    read_path : str
        Base path to ExoMol database directory
    T_list : list of float
        List of temperatures in Kelvin
    states_df : pandas.DataFrame, optional
        Complete states data used to calculate Q(T) when the requested
        temperature is absent from the partition-function file.

    Returns
    -------
    np.ndarray
        Partition function values Q(T) for each temperature, shape (n_temps,)
    """
    pf_filename = read_path + '/'.join(data_info) + '/' + '__'.join(data_info[-2:]) + '.pf'
    pf_col_name = ['T', 'Q']
    try:
        # Use Python engine for regex separator; avoids pandas warning and parses correctly.
        pf_df = pd.read_csv(pf_filename, sep=r'\s+', names=pf_col_name, header=None, engine='python')
    except:
        raise ValueError('No partition function file. Please check the file path.')
    from ..calculation.calculate_partition_func import partition_values

    def load_states():
        full_states = states_df() if callable(states_df) else states_df
        if full_states is None:
            filename = get_statesfile(read_path, data_info)
            full_states = pd.read_csv(
                filename, sep=r'\s+', header=None, usecols=[1, 2],
                names=['E', 'g'], dtype=float,
            )
        return full_states['E'].values, full_states['g'].values

    return partition_values(pf_df, T_list, load_states)
