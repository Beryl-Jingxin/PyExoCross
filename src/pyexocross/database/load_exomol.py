"""
Load ExoMol database files.

This module provides functions for reading ExoMol states files, transition files,
partition function files, and broadening parameter files.
"""
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

BROAD_DEFAULT_GAMMA_L = 0.07
BROAD_DEFAULT_N_AIR = 0.5
BROAD_RECIPE_QN_COLUMNS = {
    'a0': ('Jpp', None),
    'j': ('Jpp', None),
    'a1': ('Jpp', 'Jp'),
    'jj': ('Jpp', 'Jp'),
    'm0': ('abs_m', None),
    'm1': ('m', None),
    'k0': ('Kp', None),
    'v0': ('vp', None),
    'k1': ('Jpp', 'Kp'),
    'v1': ('Jpp', 'vp'),
}

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
    t.end()     
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
    QNcolumns = [col for col in list(dict.fromkeys(['id','E','g','J'] + optional_cols + QNslabel_list + QNs_label)) if col in available_cols]

    # LTE 
    if NLTEMethod == 'L':
        states_part_df = states_parts_df[QNcolumns]
    # Non-LTE using two temperatures
    elif NLTEMethod == 'T':
        QNcolumns_2Ts = ['id','E','g','J'] + optional_cols + QNs_label + rot_label + vib_label
        # Remove duplicates while preserving order
        QNcolumns_2T = [col for col in list(dict.fromkeys(QNcolumns_2Ts)) if col in available_cols]
        states_parts_df = states_parts_df[QNcolumns_2T]
        min_value_list = states_parts_df.sort_values('E').iloc[0]
        J_min = min_value_list['J']
        Evib_df = states_parts_df[states_parts_df[rot_label].isin(list(min_value_list[rot_label])).all(axis=1)][vib_label]
        Evib_df['Evib'] = states_parts_df['E']
        states_part_df = dd.merge(states_parts_df, Evib_df, left_on=vib_label, right_on=vib_label, how='left').dropna()
        states_part_df['Erot'] = states_part_df['E'] - states_part_df['Evib']
        states_part_df = states_part_df[QNcolumns + ['Erot', 'Evib']]
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
    folder = read_path + '/'.join(data_info) + '/'
    trans_filepaths_all = preferred_files(folder, '.trans')
    num_transfiles_all = len(trans_filepaths_all)    # The number of all transitions files including the older version files.
    trans_filepaths = []    # The list of the lastest transitions files.
    all_decompress_num = 0
    decompress_num = 0
    for i in range(num_transfiles_all):
        split_version = trans_filepaths_all[i].split('__')[-1].split('.')[0].split('_')    # Split the filenames.
        num = len(split_version)
        # There are four format filenames.
        # The lastest transitions files named in four formats:
        # 1. Filenames are named with the name of isotopologue and dataset. 
        #    End with .trans.bz2.
        #    e.g. 14N-16O__XABC.trans.bz2'
        # 2. Filenames are named with the name of isotopologue and dataset. 
        #    Also have the range of wavenumbers xxxxx-yyyyy.
        #    End with .trans.bz2.
        #    e.g. 1H2-16O__POKAZATEL__00000-00100.trans.bz2
        # 3. The older version transitions files are named with vn(version number) based on the first format of the lastest files.
        #    e.g. 14N-16O__XABC_v2.trans.bz2
        # 4. The older version transitions files are named with updated date (yyyymmdd).
        #    e.g. 1H3_p__MiZATeP__20170330.trans.bz2
        # After split the filenames:
        # The first format filenames only leave the dataset name, e.g. XABC.
        # The second format filenames only leave the range of the wavenumber, e.g. 00000-00100.
        # The third format filenames leave two parts(dataset name and version number), e.g. XABC and v2.
        # The fourth format filenames only leave the updated date, e.g. 20170330.
        # This program only process the lastest data, so extract the filenames named by the first two format.
        if num == 1:     
            original_path = trans_filepaths_all[i]
            file_size_bytes = os.path.getsize(original_path)
            trans_filepath = original_path
            if prepare and original_path.endswith('.bz2') and file_size_bytes >= LARGE_TRANS_FILE_BYTES:
                (trans_filepath, num_dec) = command_decompress(original_path)
                all_decompress_num += 1
                decompress_num += num_dec
            if split_version[0] == data_info[-1]:
                trans_filepaths.append(trans_filepath)
            elif len(split_version[0].split('-')) == 2:
                trans_filepaths.append(trans_filepath)
            else:
                pass
        else:
            pass
    print('{:45s} : {}'.format('Number of all transitions files', num_transfiles_all))
    print('{:45s} : {}'.format('Number of all decompressed transitions files', all_decompress_num))
    print('{:45s} : {}'.format('Number of new decompressed transitions files', decompress_num))
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
    folder = read_path + '/'.join(data_info) + '/'
    trans_filepaths_all = preferred_files(folder, '.trans')
    num_transfiles_all = len(trans_filepaths_all)    # The number of all transitions files including the older version files.
    trans_filepaths = []    # The list of the lastest transitions files.
    all_decompress_num = 0
    decompress_num = 0
    for i in range(num_transfiles_all):
        split_version = trans_filepaths_all[i].split('__')[-1].split('.')[0].split('_')    # Split the filenames.
        num = len(split_version)
        # There are four format filenames.
        # The lastest transitions files named in four formats:
        # 1. Filenames are named with the name of isotopologue and dataset. 
        #    End with .trans.bz2.
        #    e.g. 14N-16O__XABC.trans.bz2'
        # 2. Filenames are named with the name of isotopologue and dataset. 
        #    Also have the range of wavenumbers xxxxx-yyyyy.
        #    End with .trans.bz2.
        #    e.g. 1H2-16O__POKAZATEL__00000-00100.trans.bz2
        # 3. The older version transitions files are named with vn(version number) based on the first format of the lastest files.
        #    e.g. 14N-16O__XABC_v2.trans.bz2
        # 4. The older version transitions files are named with updated date (yyyymmdd).
        #    e.g. 1H3_p__MiZATeP__20170330.trans.bz2
        # After split the filenames:
        # The first format filenames only leave the dataset name, e.g. XABC.
        # The second format filenames only leave the range of the wavenumber, e.g. 00000-00100.
        # The third format filenames leave two parts(dataset name and version number), e.g. XABC and v2.
        # The fourth format filenames only leave the updated date, e.g. 20170330.
        # This program only process the lastest data, so extract the filenames named by the first two formats.
        if num == 1:     
            if split_version[0] == data_info[-1]:
                trans_filepath = trans_filepaths_all[i]
                file_size_bytes = os.path.getsize(trans_filepath)
                if prepare and trans_filepath.endswith('.bz2') and file_size_bytes >= LARGE_TRANS_FILE_BYTES:
                    (trans_filepath, num_dec) = command_decompress(trans_filepath)
                    all_decompress_num += 1
                    decompress_num += num_dec
                trans_filepaths.append(trans_filepath)
            elif len(split_version[0].split('-')) == 2:
                lower = int(split_version[0].split('-')[0])
                upper = int(split_version[0].split('-')[1])
                if ((lower <= int(min_wn) < upper) or 
                    (lower >= int(min_wn) and upper <= int(max_wn)) or 
                    (lower <= int(max_wn) < upper)):
                    file_size_bytes = os.path.getsize(trans_filepaths_all[i])
                    trans_filepath = trans_filepaths_all[i]
                    if prepare and trans_filepath.endswith('.bz2') and file_size_bytes >= LARGE_TRANS_FILE_BYTES:
                        (trans_filepath, num_dec) = command_decompress(trans_filepath)
                        all_decompress_num += 1
                        decompress_num += num_dec
                    trans_filepaths.append(trans_filepath)
            else:
                pass
        else:
            pass
    print('{:45s} : {}'.format('Number of all transitions files', num_transfiles_all))
    print('{:45s} : {}'.format('Number of selected transitions files', len(trans_filepaths)))
    print('{:45s} : {}'.format('Number of all decompressed transitions files', all_decompress_num))
    print('{:45s} : {}'.format('Number of new decompressed transitions files', decompress_num))
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
def read_exomol_pf(read_path, data_info, T_list):
    """
    Read ExoMol partition function file and return Q(T) for specified temperatures.

    Parameters
    ----------
    read_path : str
        Base path to ExoMol database directory
    T_list : list of float
        List of temperatures in Kelvin

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
    try:
        Q_list = pf_df[pf_df['T'].isin(T_list)]['Q']
    except:
        raise ValueError('No specified temperature dependent partition funtion value.', 
                          'Please change the temperature(s) or calculate the partition function at first.')
    if Q_list.empty:
        raise ValueError('No specified temperature dependent partition funtion value.',
                         'Please change the temperature(s) or calculate the partition function at first.')
    Q_arr = Q_list.to_numpy(dtype=float)
    return Q_arr

# Read Broadening File
def broad_recipe_qn_columns(code):
    """
    Return the physical meaning of q1/q2 for an ExoMol .broad recipe.

    The first three .broad columns are always recipe code, gamma_L, and n_air.
    The returned names describe the 4th and 5th columns, when present.
    Unknown recipes return generic q1/q2 names and are not interpreted by
    extract_broad.
    """
    recipe = str(code).lower().strip()[:2]
    return BROAD_RECIPE_QN_COLUMNS.get(recipe, ('q1', 'q2'))


def broad_required_line_columns(broad_dfs, available_columns):
    """
    Return only the line-list columns needed by the present broad recipes.

    This is deliberately narrower than keeping all quantum-number columns:
    QN filters run before line-list trimming, while broadening only needs the
    columns implied by the .broad recipe codes.
    """
    available = set(available_columns)
    required = ['J"'] if 'J"' in available else []
    recipes = set()
    for broad_df in broad_dfs:
        if 'code' not in broad_df.columns:
            recipes.add('a0')
            continue
        recipes.update(broad_df['code'].astype(str).str.lower().str.strip().str[:2])

    if recipes & {'a1', 'jj', 'm0', 'm1'} and "J'" in available:
        required.append("J'")

    if recipes & {'k0', 'k1'}:
        for candidate in ("K'", "k'"):
            if candidate in available:
                required.append(candidate)
                break

    if recipes & {'v0', 'v1'} and "v'" in available:
        required.append("v'")

    return list(dict.fromkeys(required))


def read_broad_file(fname):
    """
    Read an ExoMol .broad file while preserving recipe-dependent QN columns.
    """
    max_cols = 0
    with open(fname, 'r') as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            max_cols = max(max_cols, len(stripped.split()))
    if max_cols < 4:
        raise ValueError(f'Invalid broadening file {fname}: expected at least 4 columns.')

    broad_df = pd.read_csv(
        fname,
        sep=r'\s+',
        header=None,
        names=range(max_cols),
        comment='#',
        engine='c',
    )
    if broad_df.shape[1] < 4:
        raise ValueError(f'Invalid broadening file {fname}: expected at least 4 columns.')

    rename = {
        0: 'code',
        1: 'gamma_L',
        2: 'n_air',
        3: 'q1',
    }
    if broad_df.shape[1] > 4:
        rename[4] = 'q2'
    broad_df = broad_df.rename(columns=rename)

    for idx in range(5, broad_df.shape[1]):
        broad_df = broad_df.rename(columns={idx: f'q{idx - 2}'})

    broad_df['code'] = broad_df['code'].astype(str)
    broad_df['gamma_L'] = pd.to_numeric(broad_df['gamma_L'], errors='coerce')
    broad_df['n_air'] = pd.to_numeric(broad_df['n_air'], errors='coerce')
    for col in [c for c in broad_df.columns if c.startswith('q')]:
        broad_df[col] = pd.to_numeric(broad_df[col], errors='coerce')

    q2 = broad_df['q2'] if 'q2' in broad_df.columns else np.nan
    for alias in set(name for names in BROAD_RECIPE_QN_COLUMNS.values() for name in names if name):
        broad_df[alias] = np.nan
    for idx, code in broad_df['code'].items():
        q1_name, q2_name = broad_recipe_qn_columns(code)
        if q1_name in broad_df.columns:
            broad_df.at[idx, q1_name] = broad_df.at[idx, 'q1']
        if q2_name in broad_df.columns:
            broad_df.at[idx, q2_name] = q2.at[idx] if hasattr(q2, 'at') else np.nan
    return broad_df


def read_broad(read_path):
    """
    Read broadening parameter files for multiple broadeners.

    Reads broadening files for each specified broadener and ratio,
    or uses default values if 'DEF' is specified.

    Parameters
    ----------
    read_path : str
        Base path to ExoMol database directory

    Returns
    -------
    tuple
        A tuple containing:
        - broad : list of str, broadener names
        - ratio : list of float, mixing ratios
        - nbroad : int, number of broadeners
        - broad_dfs : list of pd.DataFrame, broadening DataFrames
    """
    from pyexocross.core import data_info, broadeners, ratios 

    broad_df = pd.DataFrame()
    broad_dfs = []
    broad = []
    ratio = []
    for i in range(len(ratios)):
        if ratios[i] != 0.0:
            if broadeners[i].upper()[0:3] == 'DEF':
                broad_df = pd.DataFrame([{
                    'code': 'Default',
                    'gamma_L': BROAD_DEFAULT_GAMMA_L,
                    'n_air': BROAD_DEFAULT_N_AIR,
                    'q1': np.nan,
                    'q2': np.nan,
                    'Jpp': np.nan,
                    'Jp': np.nan,
                }])
                broad_dfs.append(broad_df)
            else:
                broadener_name = str(broadeners[i])
                pattern_broadener = read_path + data_info[0] + '/**/*' + broadener_name + '.broad'
                if glob.glob(pattern_broadener, recursive=True) != []:
                    for fname_broadener in glob.glob(pattern_broadener, recursive=True):
                        broad_df = read_broad_file(fname_broadener)
                        broad_dfs.append(broad_df)
                else:
                    raise ValueError('The ' + broadener_name + ' boradening file does not exist.') 
            broad.append(broadeners[i])
            ratio.append(ratios[i])
    nbroad = len(broad)
    broad = list(i for i in broad if i==i)
    ratio = list(i for i in ratio if i==i)
    return(broad, ratio, nbroad, broad_dfs) 

def qn_number(values):
    return pd.to_numeric(values, errors='coerce').to_numpy(dtype=np.float64)


def m_int(value):
    value = float(value)
    rounded = round(value)
    if not np.isclose(value, rounded, atol=1e-8):
        return None
    return int(rounded)


def series_from_values(values, index):
    return pd.Series(np.asarray(values, dtype=np.float64), index=index)


def lookup_recipe_params(recipe_df, key_cols, st_key_values):
    table = recipe_df.dropna(subset=key_cols).copy()
    if table.empty:
        return None

    table = table.sort_values(key_cols).drop_duplicates(key_cols, keep='last')
    lookup = table.set_index(key_cols)[['gamma_L', 'n_air']]
    if len(key_cols) == 1:
        keys = pd.Index(st_key_values[0], name=key_cols[0])
    else:
        keys = pd.MultiIndex.from_arrays(st_key_values, names=key_cols)
    matched = lookup.reindex(keys)
    return matched['gamma_L'].to_numpy(), matched['n_air'].to_numpy()


def extract_a0_broad(recipe_df, st_df):
    table = recipe_df.dropna(subset=['q1']).copy()
    if table.empty:
        return None

    table['J_lookup'] = qn_number(table['q1'])
    table = table.sort_values('J_lookup').drop_duplicates('J_lookup', keep='last')
    j_table = table.set_index('J_lookup')[['gamma_L', 'n_air']]
    max_j = float(j_table.index.max())
    jpp = np.minimum(qn_number(st_df['J"']), max_j)
    matched = j_table.reindex(jpp)
    return matched['gamma_L'].to_numpy(), matched['n_air'].to_numpy()


def expand_m_recipe(recipe_df, signed):
    rows = []
    for row in recipe_df.itertuples(index=False):
        m_val = getattr(row, 'q1')
        if pd.isna(m_val):
            continue
        gamma_l = getattr(row, 'gamma_L')
        n_air = getattr(row, 'n_air')
        if pd.isna(gamma_l) or pd.isna(n_air):
            continue

        m_value = m_int(m_val)
        if m_value is None:
            continue
        if signed:
            if m_value < 0:
                jpp = abs(m_value)
                jp = jpp - 1
                rows.append((jpp, jp - jpp, gamma_l, n_air))
            elif m_value == 0:
                rows.append((0, 0, gamma_l, n_air))
            else:
                rows.append((m_value, 0, gamma_l, n_air))
                rows.append((m_value - 1, 1, gamma_l, n_air))
        else:
            if m_value <= 0:
                continue
            rows.append((m_value, -1, gamma_l, n_air))
            rows.append((m_value, 0, gamma_l, n_air))
            rows.append((m_value - 1, 1, gamma_l, n_air))

    if not rows:
        return pd.DataFrame(columns=['Jpp_lookup', 'dJ_lookup', 'gamma_L', 'n_air'])
    return pd.DataFrame(rows, columns=['Jpp_lookup', 'dJ_lookup', 'gamma_L', 'n_air'])


def extract_broad(broad_df, st_df):
    """
    Extract broadening parameters for transitions based on ExoMol recipes.

    Supports the ExoCross/ExoMol diet recipes a0/J, a1/JJ, m0, m1, K0/V0,
    and K1/V1. Recipes that require upper-state K/v columns are used only
    when matching columns are present in ``st_df``.

    Parameters
    ----------
    broad_df : pd.DataFrame
        Broadening DataFrame with columns: code, gamma_L, n_air, Jpp
    st_df : pd.DataFrame
        States DataFrame for transitions with 'J"' column

    Returns
    -------
    tuple of (pd.Series, pd.Series)
        A tuple containing:
        - gamma_L : pd.Series, Lorentzian HWHM coefficients
        - n_air : pd.Series, temperature exponents
    """
    if len(st_df) == 0:
        return (series_from_values([], st_df.index), series_from_values([], st_df.index))

    broad_df = broad_df.copy()
    if 'q1' not in broad_df.columns and 'Jpp' in broad_df.columns:
        broad_df['q1'] = broad_df['Jpp']
    if 'q2' not in broad_df.columns:
        broad_df['q2'] = np.nan
    if 'code' not in broad_df.columns:
        broad_df['code'] = 'a0'

    broad_df['code_norm'] = broad_df['code'].astype(str).str.lower().str.strip().str[:2]
    broad_df['gamma_L'] = pd.to_numeric(broad_df['gamma_L'], errors='coerce')
    broad_df['n_air'] = pd.to_numeric(broad_df['n_air'], errors='coerce')
    broad_df['q1'] = pd.to_numeric(broad_df['q1'], errors='coerce')
    broad_df['q2'] = pd.to_numeric(broad_df['q2'], errors='coerce')

    gamma = np.full(len(st_df), BROAD_DEFAULT_GAMMA_L, dtype=np.float64)
    n_air = np.full(len(st_df), BROAD_DEFAULT_N_AIR, dtype=np.float64)
    matched = np.zeros(len(st_df), dtype=bool)

    jpp = qn_number(st_df['J"'])
    if "J'" in st_df.columns:
        jp = qn_number(st_df["J'"])
    else:
        jp = jpp
    dj = jp - jpp

    recipe_priority = ['a0', 'j', 'a1', 'jj', 'm0', 'm1', 'k0', 'v0', 'k1', 'v1']
    for recipe in recipe_priority:
        recipe_df = broad_df[broad_df['code_norm'] == recipe].copy()
        if recipe_df.empty:
            continue

        result = None
        if recipe in ('a0', 'j'):
            result = extract_a0_broad(recipe_df, st_df)
        elif recipe in ('a1', 'jj') and "J'" in st_df.columns:
            recipe_df['Jpp_lookup'] = qn_number(recipe_df['q1'])
            recipe_df['dJ_lookup'] = qn_number(recipe_df['q2']) - recipe_df['Jpp_lookup']
            result = lookup_recipe_params(
                recipe_df,
                ['Jpp_lookup', 'dJ_lookup'],
                [jpp, dj],
            )
        elif recipe in ('m0', 'm1'):
            expanded = expand_m_recipe(recipe_df, signed=(recipe == 'm1'))
            result = lookup_recipe_params(
                expanded,
                ['Jpp_lookup', 'dJ_lookup'],
                [jpp, dj],
            )
        elif recipe in ('k0', 'v0', 'k1', 'v1'):
            qn_prefix = 'K' if recipe.startswith('k') else 'v'
            candidates = [
                f"{qn_prefix}'",
                f"{qn_prefix.upper()}'",
                f"{qn_prefix.lower()}'",
                "K'",
                "v'",
            ]
            qn_col = next((col for col in candidates if col in st_df.columns), None)
            if qn_col is None:
                continue
            qn = qn_number(st_df[qn_col])
            recipe_df['QN_lookup'] = qn_number(recipe_df['q1'] if recipe in ('k0', 'v0') else recipe_df['q2'])
            if recipe in ('k0', 'v0'):
                result = lookup_recipe_params(recipe_df, ['QN_lookup'], [qn])
            else:
                recipe_df['Jpp_lookup'] = qn_number(recipe_df['q1'])
                result = lookup_recipe_params(
                    recipe_df,
                    ['Jpp_lookup', 'QN_lookup'],
                    [jpp, qn],
                )

        if result is not None:
            gamma_values, n_values = result
            fill_mask = ~(np.isnan(gamma_values) | np.isnan(n_values))
            gamma[fill_mask] = gamma_values[fill_mask]
            n_air[fill_mask] = n_values[fill_mask]
            matched[fill_mask] = True

    if not matched.all():
        a0_df = broad_df[broad_df['code_norm'].isin(['a0', 'j'])]
        fallback = extract_a0_broad(a0_df, st_df) if not a0_df.empty else None
        if fallback is not None:
            gamma_values, n_values = fallback
            fill_mask = (~matched) & ~(np.isnan(gamma_values) | np.isnan(n_values))
            gamma[fill_mask] = gamma_values[fill_mask]
            n_air[fill_mask] = n_values[fill_mask]

    return(series_from_values(gamma, st_df.index), series_from_values(n_air, st_df.index))
