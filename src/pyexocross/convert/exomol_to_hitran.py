"""
Convert ExoMol format to HITRAN format.

This module provides functions for converting ExoMol states and transitions
to HITRAN format, including quantum number conversion and format translation.
"""
import glob
import numpy as np
import pandas as pd
import dask.dataframe as dd
from tqdm import tqdm
from functools import partial
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

from ..base import (
    Timer,
    ensure_dir,
    MAX_LARGE_FILE_WORKERS,
    print_conversion_info,
    print_file_info,
    Tref,
)
from ..base.large_file import (
    save_large_txt,
    is_large_trans_file,
    read_trans_chunks,
    process_large_chunks,
)
from ..calculation.calculate_para import cal_v, cal_uncertainty
from ..calculation.calculate_intensity import cal_abscoefs
from ..database.load_exomol import get_transfiles, read_exomol_pf, extract_broad

_USE_THREAD_POOL = True


def _executor_context(max_workers):
    """
    Use thread pool for conversion workers.

    The conversion pipeline relies on runtime globals initialised by
    ``Config.to_globals()``. With spawn-based process workers (default on
    macOS), child processes do not inherit those globals and conversion fails.
    Threads share the same interpreter state and keep behaviour stable.
    """
    return ThreadPoolExecutor(max_workers=max_workers)


def _ensure_conversion_globals():
    """
    Ensure legacy-style configuration globals are available in this module.
    Needed for worker processes spawned by ProcessPoolExecutor.
    """
    from pyexocross.core import ( 
        Conversion,
        ConversionMinFreq,
        ConversionMaxFreq,
        ConversionUnc,
        ConversionThreshold,
        GlobalQNLabel_list,
        GlobalQNFormat_list,
        LocalQNLabel_list,
        LocalQNFormat_list,
        check_uncertainty,
        check_lifetime,
        check_gfactor,
        QNslabel_list,
        read_path,
        save_path,
        data_info,
        ncputrans,
        ncpufiles,
        abundance,
        species_main_id,
        species_sub_id,
    )
    globals().update(
        dict(
            Conversion=Conversion,
            ConversionMinFreq=ConversionMinFreq,
            ConversionMaxFreq=ConversionMaxFreq,
            ConversionUnc=ConversionUnc,
            ConversionThreshold=ConversionThreshold,
            GlobalQNLabel_list=GlobalQNLabel_list,
            GlobalQNFormat_list=GlobalQNFormat_list,
            LocalQNLabel_list=LocalQNLabel_list,
            LocalQNFormat_list=LocalQNFormat_list,
            check_uncertainty=check_uncertainty,
            check_lifetime=check_lifetime,
            check_gfactor=check_gfactor,
            QNslabel_list=QNslabel_list,
            read_path=read_path,
            save_path=save_path,
            data_info=data_info,
            ncputrans=ncputrans,
            ncpufiles=ncpufiles,
            abundance=abundance,
            species_main_id=species_main_id,
            species_sub_id=species_sub_id,
        )
    )


def read_unc_states(states_df):
    """
    Read and filter states DataFrame based on uncertainty criteria.

    Filters states by uncertainty threshold if conversion is enabled and
    uncertainty filter is specified. Prepares states DataFrame for conversion
    by selecting relevant columns and setting up index.

    Parameters
    ----------
    states_df : pd.DataFrame
        Input states DataFrame with state information

    Returns
    -------
    pd.DataFrame
        Filtered and formatted states DataFrame ready for conversion
    """
    if Conversion != 0:
        if ConversionUnc != 'None':
            states_unc_df = states_df[states_df['unc'].astype(float) <= ConversionUnc]
        else:
            states_unc_df = states_df
        states_unc_df['id'] = pd.to_numeric(states_unc_df['id'])
        states_unc_df.set_index(['id'], inplace=True, drop=False)
    else:
        states_unc_df = states_df
        states_unc_df['id'] = pd.to_numeric(states_unc_df['id'])
        states_unc_df.set_index(['id'], inplace=True, drop=False)
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
    fullcolname = ['id','E','g','J'] + col_unc + col_lifetime + col_gfac + QNslabel_list
    states_unc_df = states_unc_df.iloc[:, : len(fullcolname)]
    states_unc_df.columns = fullcolname  
    colnames = ['id','E','g'] + col_unc + GlobalQNLabel_list + LocalQNLabel_list
    states_unc_df = states_unc_df[colnames] 
    states_unc_df = convert_QNValues_exomol2hitran(states_unc_df, GlobalQNLabel_list, LocalQNLabel_list)
    states_unc_df.index.name='index'
    return states_unc_df

## ExoMol to HITRAN
def convert_QNValues_exomol2hitran(states_unc_df, GlobalQNLabel_list, LocalQNLabel_list):
    """
    Convert quantum number values from ExoMol format to HITRAN format.

    Converts symmetry labels (Gtot, Gvib, Grot) from numeric codes to HITRAN
    string format, and converts taui from numeric to character format.

    Parameters
    ----------
    states_unc_df : pd.DataFrame
        States DataFrame with ExoMol quantum number values
    GlobalQNLabel_list : list of str
        List of global quantum number labels
    LocalQNLabel_list : list of str
        List of local quantum number labels

    Returns
    -------
    pd.DataFrame
        States DataFrame with quantum number values converted to HITRAN format
    """
    QNLabel_list = GlobalQNLabel_list+LocalQNLabel_list
    if 'Gtot' in QNLabel_list:
        states_unc_df["Gtot"] = (states_unc_df["Gtot"].replace('1',"A1'").replace('2',"A2'").replace('3',"E'")
                                 .replace('4','A1"').replace('5','A2"').replace('6','E"'))
    if 'Gvib' in QNLabel_list:
        states_unc_df["Gvib"] = (states_unc_df["Gvib"].replace('1',"A1'").replace('2',"A2'").replace('3',"E'")
                                 .replace('4','A1"').replace('5','A2"').replace('6','E"'))   
    if 'Grot' in QNLabel_list:
        states_unc_df["Grot"] = (states_unc_df["Grot"].replace('1',"A1'").replace('2',"A2'").replace('3',"E'")
                                 .replace('4','A1"').replace('5','A2"').replace('6','E"'))                   
    if 'taui' in QNLabel_list:
        states_unc_df["taui"] = states_unc_df["taui"].replace('0','s').replace('1','a')
    return states_unc_df

def linelist_ExoMol2HITRAN(states_unc_df,trans_part_df):
    """
    Convert ExoMol linelist to HITRAN format for a chunk of transitions.

    Processes transitions by merging with state information, calculating
    wavenumbers, and extracting necessary parameters for HITRAN format.

    Parameters
    ----------
    states_unc_df : pd.DataFrame
        States DataFrame with uncertainty filtering applied
    trans_part_df : pd.DataFrame
        Transition DataFrame chunk with columns ['u', 'l', 'A']

    Returns
    -------
    tuple
        A tuple containing:
        - exomolst_df : pd.DataFrame, states information for transitions
        - v : pd.Series, wavenumber values
        - A : np.ndarray, Einstein A coefficients
        - Epp : np.ndarray, lower state energies
        - gp : np.ndarray, upper state degeneracies
        - gpp : np.ndarray, lower state degeneracies
        - unc : np.ndarray, combined uncertainties
    """
    # Optimized: use indexed lookup instead of two merges to reduce memory
    if isinstance(states_unc_df, dd.DataFrame):
        states_unc_df = states_unc_df.compute()
    if isinstance(trans_part_df, dd.DataFrame):
        trans_part_df = trans_part_df.compute()
    
    # Set index for fast lookup (more memory efficient than merge)
    states_indexed = states_unc_df.set_index('id')
    
    # Filter trans_df to only include transitions where both states exist
    valid_mask = trans_part_df['u'].isin(states_indexed.index) & trans_part_df['l'].isin(states_indexed.index)
    trans_part_df = trans_part_df[valid_mask]
    
    if len(trans_part_df) == 0:
        return (pd.DataFrame(), pd.Series(dtype=float), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]))
    
    # Get upper and lower state data using vectorized lookup
    u_states = states_indexed.loc[trans_part_df['u']].reset_index(drop=True)
    l_states = states_indexed.loc[trans_part_df['l']].reset_index(drop=True)
    
    # Rename columns with suffixes (id column is already dropped by reset_index)
    u_states.columns = [col + "'" for col in u_states.columns]
    l_states.columns = [col + '"' for col in l_states.columns]
    
    # Combine trans and states data
    merged_df = pd.concat([
        trans_part_df[['A']].reset_index(drop=True),
        u_states,
        l_states
    ], axis=1)
    
    merged_df['v'] = cal_v(merged_df["E'"].values, merged_df['E"'].values)
    # Filter out transitions where Ep < Epp (v < 0), skip invalid line list entries
    merged_df = merged_df[merged_df['v'] >= 0]
    if len(merged_df) == 0:
        return (pd.DataFrame(), pd.Series(dtype=float), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]))
    merged_df = merged_df[merged_df['v'].between(ConversionMinFreq, ConversionMaxFreq)]
    v = merged_df['v']
    A = merged_df['A'].values
    Epp = merged_df['E"'].values
    gp = merged_df["g'"].values
    gpp = merged_df['g"'].values
    exomolst_df = merged_df.drop(columns=['A',"E'",'E"','v',"g'",'g"'])
    unc = cal_uncertainty(exomolst_df["unc'"], exomolst_df['unc"'])
    exomolst_df.drop(columns=["unc'",'unc"'], inplace=True)
    return (exomolst_df, v, A, Epp, gp, gpp, unc)

def broadener_ExoMol2HITRAN(exomolst_df):
    """
    Extract broadening parameters from ExoMol broadening files.

    Reads air and self-broadening files and extracts gamma_L and n_air
    parameters. Uses default values if broadening files are not found.

    Parameters
    ----------
    exomolst_df : pd.DataFrame
        States DataFrame for transitions, used to determine number of rows

    Returns
    -------
    tuple of (np.ndarray, np.ndarray, np.ndarray)
        A tuple containing:
        - gamma_air : np.ndarray, air-broadening coefficients
        - gamma_self : np.ndarray, self-broadening coefficients
        - n_air : np.ndarray, temperature exponents
    """
    broad_col_name = ['code', 'gamma_L', 'n_air', 'Jpp']
    default_broad_df = pd.DataFrame(columns=broad_col_name)
    default_gamma_L = 0.07
    default_n_air = 0.5
    default_broad_df = pd.DataFrame([['code', default_gamma_L, default_n_air,'Jpp']],columns=broad_col_name)
    air_broad_df = pd.DataFrame(columns=broad_col_name)
    rows = len(exomolst_df)
    pattern_air = read_path + data_info[0] + '/**/*air.broad'
    if glob.glob(pattern_air, recursive=True) != []:
        for fname_air in glob.glob(pattern_air, recursive=True):
            air_broad_df = pd.read_csv(fname_air, sep='\s+', names=broad_col_name, header=None, engine='python')
            gamma_air = extract_broad(air_broad_df,exomolst_df)[0].values
            n_air = extract_broad(air_broad_df,exomolst_df)[1].values
    else:
        gamma_air= np.full((1,rows),default_broad_df['gamma_L'][0])[0]
        n_air = np.full((1,rows),default_broad_df['n_air'][0])[0]
    pattern_self = read_path + data_info[0] + '/**/*self.broad'
    if glob.glob(pattern_self, recursive=True) != []:
        for fname_self in glob.glob(pattern_self, recursive=True):
            self_broad_df = pd.read_csv(fname_self, sep='\s+', names=broad_col_name, header=None, engine='python')
            gamma_self = extract_broad(self_broad_df,exomolst_df)[0].values
    else:
        gamma_self= np.full((1,rows),default_broad_df['gamma_L'][0])[0]  
    return (gamma_air, gamma_self, n_air)

def convert_QNFormat_exomol2hitran(exomolst_df, GlobalQNLabel_list, GlobalQNFormat_list, 
                                   LocalQNLabel_list, LocalQNFormat_list):
    """
    Convert quantum number formats from ExoMol to HITRAN format.

    Formats quantum number values according to HITRAN specifications,
    combining global and local quantum numbers into V'V"Q'Q" columns.

    Parameters
    ----------
    exomolst_df : pd.DataFrame
        States DataFrame with quantum number columns
    GlobalQNLabel_list : list of str
        List of global quantum number labels
    GlobalQNFormat_list : list of str
        List of format strings for global quantum numbers
    LocalQNLabel_list : list of str
        List of local quantum number labels
    LocalQNFormat_list : list of str
        List of format strings for local quantum numbers

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ["V'", 'V"', "Q'", 'Q"'] containing formatted quantum numbers
    """
    n_gQN = len(GlobalQNLabel_list)
    for i in range(n_gQN):
        gQN_format = GlobalQNFormat_list[i].replace("%",'{: >')+'}'
        gQN_label = GlobalQNLabel_list[i]
        try:
            if 'd' in gQN_format or 'f' in gQN_format: 
                exomolst_df[gQN_label+"'"] = pd.Series(pd.to_numeric(exomolst_df[gQN_label+"'"].values)).parallel_map(gQN_format.format)
                exomolst_df[gQN_label+'"'] = pd.Series(pd.to_numeric(exomolst_df[gQN_label+'"'].values)).parallel_map(gQN_format.format)
            elif 's' in gQN_format or 'a' in gQN_format: 
                exomolst_df[gQN_label+"'"] = pd.Series(exomolst_df[gQN_label+"'"].str.replace('(','').str.replace(')','')).parallel_map(gQN_format.format)
                exomolst_df[gQN_label+'"'] = pd.Series(exomolst_df[gQN_label+'"'].str.replace('(','').str.replace(')','')).parallel_map(gQN_format.format)
        except:
            if 'd' in gQN_format or 'f' in gQN_format: 
                exomolst_df[gQN_label+"'"] = pd.Series(pd.to_numeric(exomolst_df[gQN_label+"'"].values)).map(gQN_format.format)
                exomolst_df[gQN_label+'"'] = pd.Series(pd.to_numeric(exomolst_df[gQN_label+'"'].values)).map(gQN_format.format)
            elif 's' in gQN_format or 'a' in gQN_format: 
                exomolst_df[gQN_label+"'"] = pd.Series(exomolst_df[gQN_label+"'"].str.replace('(','').str.replace(')','')).map(gQN_format.format)
                exomolst_df[gQN_label+'"'] = pd.Series(exomolst_df[gQN_label+'"'].str.replace('(','').str.replace(')','')).map(gQN_format.format)
    try:
        globalQNp = exomolst_df[[GlobalQNLabel_list[i]+"'" for i in range(n_gQN)]].astype(str).sum(axis=1).parallel_map('{: >15}'.format)  
        globalQNpp = exomolst_df[[GlobalQNLabel_list[i]+'"' for i in range(n_gQN)]].astype(str).sum(axis=1).parallel_map('{: >15}'.format)    
    except:
        globalQNp = exomolst_df[[GlobalQNLabel_list[i]+"'" for i in range(n_gQN)]].astype(str).sum(axis=1).map('{: >15}'.format)  
        globalQNpp = exomolst_df[[GlobalQNLabel_list[i]+'"' for i in range(n_gQN)]].astype(str).sum(axis=1).map('{: >15}'.format)     
                
    n_lQN = len(LocalQNLabel_list)
    for i in range(n_lQN):
        lQN_format = LocalQNFormat_list[i].replace("%",'{: >')+'}'
        lQN_label = LocalQNLabel_list[i]
        try:
            if 'd' in lQN_format or 'f' in lQN_format: 
                exomolst_df[lQN_label+"'"] = pd.Series(pd.to_numeric(exomolst_df[lQN_label+"'"].values)).parallel_map(lQN_format.format)
                exomolst_df[lQN_label+'"'] = pd.Series(pd.to_numeric(exomolst_df[lQN_label+'"'].values)).parallel_map(lQN_format.format)
            elif 's' in lQN_format or 'a' in lQN_format: 
                exomolst_df[lQN_label+"'"] = pd.Series(exomolst_df[lQN_label+"'"].str.replace('(','').str.replace(')','')).parallel_map(lQN_format.format)
                exomolst_df[lQN_label+'"'] = pd.Series(exomolst_df[lQN_label+'"'].str.replace('(','').str.replace(')','')).parallel_map(lQN_format.format)
        except:
            if 'd' in lQN_format or 'f' in lQN_format: 
                exomolst_df[lQN_label+"'"] = pd.Series(pd.to_numeric(exomolst_df[lQN_label+"'"].values)).map(lQN_format.format)
                exomolst_df[lQN_label+'"'] = pd.Series(pd.to_numeric(exomolst_df[lQN_label+'"'].values)).map(lQN_format.format)
            elif 's' in lQN_format or 'a' in gQN_format: 
                exomolst_df[lQN_label+"'"] = pd.Series(exomolst_df[lQN_label+"'"].str.replace('(','').str.replace(')','')).map(lQN_format.format)
                exomolst_df[lQN_label+'"'] = pd.Series(exomolst_df[lQN_label+'"'].str.replace('(','').str.replace(')','')).map(lQN_format.format)
    try:
        localQNp = exomolst_df[[LocalQNLabel_list[i]+"'" for i in range(n_lQN)]].astype(str).sum(axis=1).parallel_map('{: >15}'.format)  
        localQNpp = exomolst_df[[LocalQNLabel_list[i]+'"' for i in range(n_lQN)]].astype(str).sum(axis=1).parallel_map('{: >15}'.format)            
    except:
        localQNp = exomolst_df[[LocalQNLabel_list[i]+"'" for i in range(n_lQN)]].astype(str).sum(axis=1).map('{: >15}'.format)  
        localQNpp = exomolst_df[[LocalQNLabel_list[i]+'"' for i in range(n_lQN)]].astype(str).sum(axis=1).map('{: >15}'.format)    
    QN_df = pd.concat([globalQNp,globalQNpp,localQNp,localQNpp],axis='columns')
    QN_df.columns = ["V'", 'V"', "Q'", 'Q"'] 
    return QN_df

def error_code(unc):
    """
    Convert uncertainty values to HITRAN error codes.

    Maps uncertainty values to 6-digit error codes according to HITRAN
    convention based on uncertainty magnitude ranges.

    Parameters
    ----------
    unc : np.ndarray
        Uncertainty values array, shape (n_levels,)

    Returns
    -------
    np.ndarray
        Error code array as integers, shape (n_levels,)
    """
    unc[(1<=unc)] = '000000'
    unc[(0.1<=unc) & (unc<1)] = '100000'
    unc[(0.01<=unc) & (unc<0.1)] = '200000'
    unc[(0.001<=unc) & (unc<0.01)] = '300000'
    unc[(0.0001<=unc) & (unc<0.001)] = '400000'
    unc[(0.00001<=unc) & (unc<0.0001)] = '500000'
    unc[(0.000001<=unc) & (unc<0.00001)] = '600000'
    unc[(0.0000001<=unc) & (unc<0.000001)] = '700000'
    unc[(0.00000001<=unc) & (unc<0.0000001)] = '800000'
    unc[(unc<0.00000001)] = '900000'
    unc = unc.astype(int)
    return unc  

def convert_exomol2hitran_linelist(states_df, trans_part_df):
    """
    Convert ExoMol linelist to HITRAN format for a transition chunk.

    Orchestrates the conversion process by calling helper functions to
    process states, transitions, broadening parameters, and quantum numbers.

    Parameters
    ----------
    states_df : pd.DataFrame
        Complete states DataFrame
    trans_part_df : pd.DataFrame
        Transition DataFrame chunk

    Returns
    -------
    pd.DataFrame
        HITRAN format linelist DataFrame
    """
    states_unc_df = read_unc_states(states_df)
    (exomolst_df, v, A, Epp, gp, gpp, unc_values) = linelist_ExoMol2HITRAN(states_unc_df,trans_part_df)
    (gamma_air, gamma_self, n_air) = broadener_ExoMol2HITRAN(exomolst_df)
    QN_df =  convert_QNFormat_exomol2hitran(exomolst_df, GlobalQNLabel_list, GlobalQNFormat_list, 
                                            LocalQNLabel_list, LocalQNFormat_list)
    Q = read_exomol_pf(read_path, data_info, [Tref])
    I = cal_abscoefs([Tref], Q, Epp, gp, A, v, abundance)[0]
    unc = error_code(unc_values) 
    nrows = len(A)
    delta_air = ['']*nrows 
    iref = ['']*nrows 
    flag = ['']*nrows 
    '''
    hitran_column_name = ['M','I','v','S','A','gamma_air','gamma_self',
                          'E"','n_air','delta_air','Vp','Vpp','Qp','Qpp',
                          'Ierr','Iref','flag','gp','gpp']
    '''
    hitran_begin_dic = {
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
    }
    hitran_begin_df = pd.DataFrame(hitran_begin_dic)
    hitran_end_dic = {'Ierr':unc,'Iref':iref,'flag':flag,"g'":gp, 'g"':gpp}
    hitran_end_df = pd.DataFrame(hitran_end_dic)

    hitran_res_df = pd.concat([hitran_begin_df, QN_df, hitran_end_df], axis='columns')
    if ConversionThreshold != 'None':
        hitran_res_df = hitran_res_df[hitran_res_df['S'] >= ConversionThreshold]
    # hitran_res_df = hitran_res_df.sort_values('v')
    return hitran_res_df

def process_exomol2hitran_linelist(states_df, trans_filepath):
    """
    Process a single ExoMol transition file and convert to HITRAN format.

    Handles both large and small files, using appropriate processing strategies
    for memory efficiency.

    Parameters
    ----------
    states_df : pd.DataFrame
        Complete states DataFrame
    trans_filepath : str
        Path to the ExoMol transition file to process

    Returns
    -------
    pd.DataFrame
        HITRAN format linelist DataFrame for this file
    """
    global _USE_THREAD_POOL
    _ensure_conversion_globals()
    trans_filename = trans_filepath.split('/')[-1]
    print('Processeing transitions file:', trans_filename)
    use_cols = [0,1,2]
    use_names = ['u','l','A']
    large_file = is_large_trans_file(trans_filepath)
    trans_reader = read_trans_chunks(trans_filepath, use_cols, use_names)
    desc = 'Processing ' + trans_filename + (' (limited streaming)' if large_file else '')
    zero_factory = lambda: pd.DataFrame()
    def combine_fn(results):
        if not results:
            return zero_factory()
        return pd.concat(results, ignore_index=True)
    handler = partial(convert_exomol2hitran_linelist, states_df)
    if large_file:
        print('Large transition file detected (>1 GB). Using bounded parallel streaming to reduce memory usage.')
        try:
            hitran_res_df = process_large_chunks(
                trans_reader,
                handler,
                combine_fn,
                zero_factory,
                desc,
                max_workers=max(1, min(ncputrans, MAX_LARGE_FILE_WORKERS)),
                executor_class=ThreadPoolExecutor if _USE_THREAD_POOL else ProcessPoolExecutor,
            )
        except PermissionError:
            # Fallback for restricted environments (e.g. sandboxed multiprocessing)
            _USE_THREAD_POOL = True
            hitran_res_df = process_large_chunks(
                trans_reader,
                handler,
                combine_fn,
                zero_factory,
                desc,
                max_workers=max(1, min(ncputrans, MAX_LARGE_FILE_WORKERS)),
                executor_class=ThreadPoolExecutor,
            )
    else:
        trans_chunks = list(trans_reader)
        if len(trans_chunks) == 0:
            hitran_res_df = zero_factory()
        else:
            with _executor_context(max_workers=ncputrans) as trans_executor:
                futures = [trans_executor.submit(convert_exomol2hitran_linelist, states_df, chunk)
                           for chunk in tqdm(trans_chunks, desc=desc)]
                hitran_res_df = combine_fn([future.result() for future in tqdm(futures, desc='Combining '+trans_filename)])
    return hitran_res_df

def conversion_exomol2hitran(states_df):
    """
    Main function to convert ExoMol database format to HITRAN format.

    Processes all transition files, converts states and transitions data,
    and saves the result in HITRAN .par file format.

    Parameters
    ----------
    states_df : pd.DataFrame
        Complete states DataFrame from ExoMol database
    """
    _ensure_conversion_globals()
    print('Convert data format from ExoMol to HITRAN.')  
    print_conversion_info(ConversionMinFreq, ConversionMaxFreq, GlobalQNLabel_list, GlobalQNFormat_list, 
                          LocalQNLabel_list, LocalQNFormat_list, ConversionUnc, ConversionThreshold)
    tot = Timer()
    tot.start()
    print('Reading transitions and converting data format from ExoMol to HITRAN ...')    
    trans_filepaths = get_transfiles(read_path, data_info)
    # Process multiple files in parallel
    with _executor_context(max_workers=ncpufiles) as executor:
        # Submit reading tasks for each file
        futures = [executor.submit(process_exomol2hitran_linelist,states_df,
                                   trans_filepath) for trans_filepath in tqdm(trans_filepaths)]
        hitran_res_df = pd.concat([future.result() for future in futures])
    hitran_res_df.sort_values(by=['v'], ascending=True, inplace=True)  
    tot.end()
    print('Finished reading all transitions and converting data format from ExoMol to HITRAN!\n')

    print('Saving HITRAN format data into file ...')   
    ts = Timer()    
    ts.start()
    conversion_folder = save_path + 'conversion/ExoMol2HITRAN/'
    ensure_dir(conversion_folder) 
    conversion_path = conversion_folder + '__'.join(data_info[-2:]) + '.par'
    hitran_format = "%2s%1s%12.6f%10.3E%10.3E%5.3f%5.3f%10.4f%4.2f%8s%15s%15s%15s%15s%6s%12s%1s%7.1f%7.1f"
    save_large_txt(conversion_path, hitran_res_df, fmt=hitran_format)
    ts.end()
    hitran_fmt_list = ['%'+i for i in hitran_format.split('%')][1:]
    print_file_info('Converted HITRAN par', hitran_res_df.columns, hitran_fmt_list)
    print('Converted HITRAN par file has been saved:', conversion_path)
    print('Converted HITRAN par file has been saved!\n')  
    print('Finished converting data format from ExoMol to HITRAN!\n')
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')
