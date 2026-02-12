"""
Save ExoMol lifetimes.

This module provides functions for calculating and saving radiative lifetimes.
"""
import re
import bz2
import numpy as np
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from src.base.utils import Timer, ensure_dir
from src.base.log import print_file_info
from src.base.large_file import read_trans_chunks
from src.calculation.calculate_lifetime import cal_lifetime
from src.database.load_exomol import get_transfiles

_USE_THREAD_POOL = False


def _executor_context(max_workers):
    """
    Prefer process pool, fall back to threads if unavailable.
    """
    global _USE_THREAD_POOL
    if _USE_THREAD_POOL:
        return ThreadPoolExecutor(max_workers=max_workers)
    try:
        return ProcessPoolExecutor(max_workers=max_workers)
    except PermissionError:
        _USE_THREAD_POOL = True
        return ThreadPoolExecutor(max_workers=max_workers)

# Process Lifetime
def process_exomol_lifetime(states_df, trans_filepath):
    """
    Process a single transition file to calculate radiative lifetimes.

    Calculates lifetime for each state from transition probabilities.

    Parameters
    ----------
    states_df : pd.DataFrame
        States DataFrame with 'id' column
    trans_filepath : str
        Path to the transition file to process

    Returns
    -------
    np.ndarray
        Lifetime array corresponding to states in states_df, shape (n_states,)
    """
    # Ensure legacy-style globals are available in worker processes.
    from pyexocross.core import ncputrans  

    global _USE_THREAD_POOL
    trans_filename = trans_filepath.split('/')[-1]
    print('Processeing transitions file:', trans_filename)
    use_cols = [0,1,2]
    use_names = ['u','l','A']
    trans_reader = read_trans_chunks(trans_filepath, use_cols, use_names)
    desc = 'Processing ' + trans_filename
    trans_chunks = list(trans_reader)
    if len(trans_chunks) == 0:
        lifetime = np.zeros(len(states_df))
    else:
        with _executor_context(max_workers=ncputrans) as trans_executor:
            futures = [trans_executor.submit(cal_lifetime, states_df, chunk)
                        for chunk in tqdm(trans_chunks, desc=desc)]
            lifetime = np.sum([future.result() for future in tqdm(futures, desc='Combining '+trans_filename)], axis=0)
    return lifetime

def insert_exomol_lifetime_column(states_df, lifetime, states_col, states_fmt):
    """
    Insert lifetime column into states DataFrame at appropriate position.

    Parameters
    ----------
    states_df : pd.DataFrame
        States DataFrame
    lifetime : np.ndarray
        Lifetime array, shape (n_states,)
    states_col : list of str
        Column names for states file
    states_fmt : list of str
        Format strings for each column

    Returns
    -------
    tuple of (pd.DataFrame, list) 
        A tuple containing:
        - states_lifetime_df : pd.DataFrame, states DataFrame with lifetime column
        - states_lifetime_fmt : list of str, updated format strings
    """
    # Import legacy-style configuration variables from core (set via Config.to_globals()).
    from pyexocross.core import check_uncertainty, check_lifetime  # type: ignore

    states_lifetime_df = states_df.copy()
    states_lifetime_df.columns = states_col
    if check_uncertainty == True:
        insert_at = 5      # Default: after 'Unc' (id, E, g, J, unc, ...)
    else:
        insert_at = 4      # Default: after 'J' (id, E, g, J, ...)
    if check_lifetime == True or states_lifetime_df.columns.str.contains('tau').any() == True:
        try:
            states_lifetime_df = states_lifetime_df.drop(columns=['tau'], axis=1)
        except:
            states_lifetime_df = states_lifetime_df.drop(columns=[insert_at], axis=1)
        states_lifetime_fmt = states_fmt[:insert_at] + ['%12s'] + states_fmt[insert_at+1:] 
    else:
        states_lifetime_fmt = states_fmt[:insert_at] + ['%12s'] + states_fmt[insert_at:]   
    states_lifetime_df.insert(insert_at, 'tau', lifetime)     
    return (states_lifetime_df, states_lifetime_fmt)

def convert_dtype_by_format(df, fmt_list):
    """
    Convert DataFrame column dtypes based on format strings.

    Maps format strings (%d, %f, %s) to appropriate pandas dtypes.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to convert
    fmt_list : list of str
        Format strings for each column (e.g., ['%12d', '%12.6f', '%15s'])

    Returns
    -------
    pd.DataFrame
        DataFrame with converted dtypes
    """
    for i, fmt in enumerate(fmt_list):
        if i < len(df.columns):
            col = df.columns[i]
            if re.search(r'%[\d.]*d', fmt):  
                df[col] = pd.to_numeric(df[col], errors='coerce').astype(int)
            elif re.search(r'%[\d.]*[fFeEgG]', fmt):  # Match f, F, e, E, g, G (all float formats)
                df[col] = pd.to_numeric(df[col], errors='coerce').astype(float)
            elif re.search(r'%[\d.]*s', fmt):  
                df[col] = df[col].astype('string')
    return df

def save_exomol_lifetime(read_path, states_df, states_col, states_fmt):
    """
    Main function to calculate and save radiative lifetimes for ExoMol database.

    Processes all transition files, calculates lifetimes, inserts lifetime column
    into states DataFrame, and saves results in .states format file.

    Parameters
    ----------
    read_path : str
        Base path to ExoMol database directory
    states_df : pd.DataFrame
        Complete states DataFrame
    states_col : list of str
        Column names for states file
    states_fmt : list of str
        Format strings for each column
    """
    np.seterr(divide='ignore', invalid='ignore')
    print('Calculate lifetimes.')  
    # Import legacy-style configuration variables from core (set via Config.to_globals()).
    from pyexocross.core import data_info, ncpufiles, save_path, CompressYN  # type: ignore

    t = Timer()
    t.start()
    print('Reading all transitions and calculating lifetimes ...')
    trans_filepaths = get_transfiles(read_path, data_info)
    # Process multiple files in parallel
    with _executor_context(max_workers=ncpufiles) as executor:
        # Submit reading tasks for each file
        futures = [executor.submit(process_exomol_lifetime, states_df, 
                                   trans_filepath) for trans_filepath in tqdm(trans_filepaths)]
        lifetime_result = 1 / sum(np.array([future.result() for future in futures]))
        lifetime = np.array([f'{x:>12.4E}'.replace('INF','Inf') for x in lifetime_result])
    t.end()
    print('Finished reading all transitions and calculating lifetimes!\n')

    print('Saving lifetimes into file ...')   
    ts = Timer()    
    ts.start()  
    
    # Insert lifetime column into states_df
    (states_lifetime_df, states_lifetime_fmt) = insert_exomol_lifetime_column(states_df, lifetime, states_col, states_fmt)
    states_lifetime_df = convert_dtype_by_format(states_lifetime_df, states_lifetime_fmt)

    lf_folder = save_path + 'lifetime/'
    ensure_dir(lf_folder)

    if CompressYN == 'Y':
        lf_path = lf_folder + '__'.join(data_info[-2:]) + '.states.bz2'
        with bz2.open(lf_path, 'wt') as f:
            np.savetxt(f, states_lifetime_df, fmt=' '.join(states_lifetime_fmt))
    else:
        lf_path = lf_folder + '__'.join(data_info[-2:]) + '.states'
        np.savetxt(lf_path, states_lifetime_df, fmt=' '.join(states_lifetime_fmt))
        
    ts.end()
    states_lifetime_cols = list(states_lifetime_df.columns)
    tau_idx = states_lifetime_cols.index('tau') 
    states_lifetime_cols = states_lifetime_cols[:tau_idx] + ['tau'] + states_lifetime_cols[tau_idx:]
    states_lifetime_fmt[tau_idx] = '%12.4E'
    print_file_info('Lifetimes', list(dict.fromkeys(states_lifetime_cols)), states_lifetime_fmt)
    print('Lifetimes file has been saved:', lf_path, '\n') 
    print('Lifetimes have been saved!\n')    
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')
