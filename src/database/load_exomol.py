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

# Read Input Files
## Read ExoMol Database Files
### Read States File
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
    states_df = pd.DataFrame()
    states_filename = read_path + '/'.join(data_info) + '/' + '__'.join(data_info[-2:]) + '.states.bz2'
    if os.path.exists(states_filename):    
        chunks = pd.read_csv(states_filename, compression='bz2', sep='\s+', header=None,
                             chunksize=100_000, iterator=True, low_memory=False, dtype=object)
    elif os.path.exists(states_filename.replace('.bz2','')):
        chunks = pd.read_csv(states_filename.replace('.bz2',''), sep='\s+', header=None,
                             chunksize=100_000, iterator=True, low_memory=False, dtype=object)
    else:
        raise ValueError("No such states file, please check the read path and states filename format!")

    for chunk in chunks:
        states_df = pd.concat([states_df, chunk])
    if check_uncertainty == True:
        states_df = states_df.rename(columns={0:'id',1:'E',2:'g',3:'J',4:'unc'})
        # Use float64 for Unc column to avoid precision issues with very small values
        convert_dict = {'id':np.int32,'E':np.float64,'g':np.int32,'J':np.float16,'unc':np.float64}
    else:      
        states_df = states_df.rename(columns={0:'id',1:'E',2:'g',3:'J'})  
        convert_dict = {'id':np.int32,'E':np.float64,'g':np.int32,'J':np.float16}
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
    if UncFilter != 'None' :
        if check_uncertainty == True:
            states_parts_df = states_df[states_df['unc'] <= UncFilter]
            # states_parts_df.set_index(['id'], inplace=True, drop=False)
            # states_parts_df['id'] = pd.to_numeric(states_parts_df['id'])
        else:
            raise ValueError("No uncertainties in states file. Please do not use uncertainty filter.")  
    else:
        states_parts_df = states_df
        # states_parts_df.set_index(['id'], inplace=True, drop=False)
        # states_parts_df['id'] = pd.to_numeric(states_parts_df['id'])
    if check_uncertainty == True:
        col_unc = ['unc']
    else:
        col_unc = []
    if check_lifetime == True:
        col_lifetime = ['tau']
    else:
        col_lifetime = []
    if check_gfactor == True:
        col_gfac = ['gfac']
    else:
        col_gfac = []
    colname = ['id','E','g','J'] + col_unc + col_lifetime + col_gfac + QNslabel_list
    # Keep only the required columns and avoid SettingWithCopyWarning by working on a copy
    states_parts_df = states_parts_df.loc[:, states_parts_df.columns[:len(colname)]].copy()
    states_parts_df.columns = colname
    QNcolumns = ['id','E','g','J'] + col_unc + col_lifetime + col_gfac + QNs_label

    # LTE 
    if NLTEMethod == 'L':
        states_part_df = states_parts_df[QNcolumns]
    # Non-LTE using two temperatures
    elif NLTEMethod == 'T':
        QNcolumns_2Ts = ['id','E','g','J'] + col_unc + col_lifetime + col_gfac + QNs_label + rot_label + vib_label
        # Remove duplicates while preserving order
        QNcolumns_2T = list(dict.fromkeys(QNcolumns_2Ts))
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
    states_part_df.index.name='index'
    try:
        states_part_df['tau'] = states_part_df['tau'].astype('float')
        states_part_df['gfac'] = states_part_df['gfac'].astype('float')
    except:
        pass
    return states_part_df

### Get transitions File
def get_transfiles(read_path, data_info):
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
    trans_filepaths_all = glob.glob((read_path + '/'.join(data_info) + '/' + '*.trans.bz2'))
    if trans_filepaths_all == []:
        trans_filepaths_all = glob.glob((read_path + '/'.join(data_info) + '/' + '*.trans'))
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
            if original_path.endswith('.bz2') and file_size_bytes >= LARGE_TRANS_FILE_BYTES:
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

def get_part_transfiles(read_path, data_info):
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

    Returns
    -------
    list of str
        List of transition file paths within the specified wavenumber range
    """
    # Get all the transitions files from the folder including the older version files which are named by vn(version number).
    trans_filepaths_all = glob.glob(read_path + '/'.join(data_info) + '/' + '*.trans.bz2')
    if trans_filepaths_all == []:
        trans_filepaths_all = glob.glob(read_path + '/'.join(data_info) + '/' + '*.trans')
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
                trans_filepaths.append(trans_filepaths_all[i])
            elif len(split_version[0].split('-')) == 2:
                lower = int(split_version[0].split('-')[0])
                upper = int(split_version[0].split('-')[1])
                if ((lower <= int(min_wn) < upper) or 
                    (lower >= int(min_wn) and upper <= int(max_wn)) or 
                    (lower <= int(max_wn) < upper)):
                    file_size_bytes = os.path.getsize(trans_filepaths_all[i])
                    trans_filepath = trans_filepaths_all[i]
                    if trans_filepath.endswith('.bz2') and file_size_bytes >= LARGE_TRANS_FILE_BYTES:
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
        pf_df = pd.read_csv(pf_filename, sep='\\s+', names=pf_col_name, header=None)
    except:
        raise ValueError('No partition function file. Please check the file path.')
    try:
        Q_list = pf_df[pf_df['T'].isin(T_list)]['Q']
    except:
        raise ValueError('No specified temperature dependent partition funtion value.', 
                          'Please change the temperature(s) or calculate the partition function at first.')
    Q_arr = Q_list.to_numpy(dtype=float)
    return Q_arr

# Read Broadening File
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
                default_gamma_L = 0.07
                default_n_air = 0.5
                broad_df = pd.DataFrame([['code', default_gamma_L, default_n_air,'Jpp']])
                broad_df = broad_df.rename(columns={0:'code', 1:'gamma_L', 2:'n_air', 3:'Jpp'})
                broad_dfs.append(broad_df)
            else:
                broadener_name = str(broadeners[i])
                pattern_broadener = read_path + data_info[0] + '/**/*' + broadener_name + '.broad'
                if glob.glob(pattern_broadener, recursive=True) != []:
                    for fname_broadener in glob.glob(pattern_broadener, recursive=True):
                        broad_df = pd.read_csv(fname_broadener, sep='\s+', header=None, engine='python')
                        broad_df = broad_df.rename(columns={0:'code', 1:'gamma_L', 2:'n_air', 3:'Jpp'})
                        broad_dfs.append(broad_df)
                else:
                    raise ValueError('The ' + broadener_name + ' boradening file does not exist.') 
            broad.append(broadeners[i])
            ratio.append(ratios[i])
    nbroad = len(broad)
    broad = list(i for i in broad if i==i)
    ratio = list(i for i in ratio if i==i)
    return(broad, ratio, nbroad, broad_dfs) 

def extract_broad(broad_df, st_df):
    """
    Extract broadening parameters for transitions based on lower state J.

    Maps broadening coefficients from broadening DataFrame to transitions
    based on lower state J quantum number.

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
    max_broad_J = max(broad_df['Jpp'])
    Jpp = st_df['J"'].values
    Jpp[Jpp > max_broad_J] = max_broad_J
    id_broad = (st_df['J"']-0.1).round(0).astype('int32').values
    gamma_L = broad_df['gamma_L'][id_broad]
    n_air = broad_df['n_air'][id_broad]
    return(gamma_L, n_air)

