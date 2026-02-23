"""
Save ExoMol cooling functions.

This module provides functions for calculating and saving cooling functions.
"""
import os
import numpy as np
import pandas as pd
import dask.dataframe as dd
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pyexocross.base.utils import Timer, ensure_dir
from pyexocross.base.log import log_tqdm, print_file_info
from pyexocross.base.large_file import (
    is_large_trans_file,
    read_trans_chunks,
)
from pyexocross.calculation.calculate_para import cal_v
from pyexocross.calculation.calculate_cooling_func import cal_cooling_func
from pyexocross.database.load_exomol import get_transfiles, read_exomol_pf

_USE_THREAD_POOL = True


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

# Process Cooling Function
def process_exomol_cooling_func_chunk(states_df,Ts,trans_df):
    """
    Process a chunk of transitions to calculate cooling functions.

    Calculates cooling function contribution from transitions at multiple temperatures.

    Parameters
    ----------
    states_df : pd.DataFrame
        States DataFrame with 'id', 'E', 'g' columns
    Ts : np.ndarray
        Temperature array
    trans_df : pd.DataFrame
        Transition DataFrame chunk with columns ['u', 'l', 'A']

    Returns
    -------
    np.ndarray
        Cooling function array, shape (n_temps,)
    """
    # Import legacy-style configuration variables from core (set via Config.to_globals()).
    from pyexocross.core import read_path, data_info, Tmax 

    # Optimized: use indexed lookup instead of two merges to reduce memory
    if isinstance(states_df, dd.DataFrame):
        states_df = states_df.compute()
    if isinstance(trans_df, dd.DataFrame):
        trans_df = trans_df.compute()
    
    # Set index for fast lookup (more memory efficient than merge)
    states_indexed = states_df.set_index('id')
    
    # Filter trans_df to only include transitions where both states exist
    valid_mask = trans_df['u'].isin(states_indexed.index) & trans_df['l'].isin(states_indexed.index)
    trans_df = trans_df[valid_mask]
    
    if len(trans_df) == 0:
        return np.zeros(len(Ts))
    
    # Use vectorized lookup instead of merge (much faster and less memory)
    Ep = states_indexed.loc[trans_df['u'], 'E'].values
    Epp = states_indexed.loc[trans_df['l'], 'E'].values
    gp = states_indexed.loc[trans_df['u'], 'g'].values
    A = trans_df['A'].values
    
    v = cal_v(Ep, Epp)
    # Filter out transitions where Ep < Epp (v < 0), skip invalid line list entries
    valid_v_mask = v >= 0
    if not valid_v_mask.any():
        return np.zeros(len(Ts))
    
    # Apply filter to all arrays
    Ep = Ep[valid_v_mask]
    Epp = Epp[valid_v_mask]
    gp = gp[valid_v_mask]
    A = A[valid_v_mask]
    v = v[valid_v_mask]
    Qs = read_exomol_pf(read_path, data_info, Ts)
    
    num = len(v)
    if num > 0:
        cooling_func = [cal_cooling_func(A, v, Ep, gp, Ts[i], Qs[i]) 
                        for i in log_tqdm(range(Tmax), desc='Calculating')]
    else:
        cooling_func = np.zeros(len(Ts))
    return cooling_func

def process_exomol_cooling_func(states_df, Ts, trans_filepath):
    """
    Process a single transition file to calculate cooling functions.

    Handles both large and small files, using appropriate processing strategies.

    Parameters
    ----------
    states_df : pd.DataFrame
        States DataFrame with 'id', 'E', 'g' columns
    Ts : np.ndarray
        Temperature array
    trans_filepath : str
        Path to the transition file to process

    Returns
    -------
    np.ndarray
        Combined cooling function array from all chunks, shape (n_temps,)
    """
    # Import legacy-style configuration variables from core (set via Config.to_globals()).
    from pyexocross.core import ncputrans 

    trans_filename = trans_filepath.split('/')[-1]
    print('Processeing transitions file:', trans_filename)
    use_cols = [0,1,2]
    use_names = ['u','l','A']
    large_file = is_large_trans_file(trans_filepath)
    trans_reader = read_trans_chunks(trans_filepath, use_cols, use_names)
    desc = 'Processing ' + trans_filename + (' (streaming)' if large_file else '')
    if large_file:
        print('Large transition file detected (>1 GB). Streaming chunks sequentially to reduce memory usage.')
        cooling_accum = None
        for trans_df_chunk in log_tqdm(trans_reader, desc=desc):
            chunk_cf = process_exomol_cooling_func_chunk(states_df, Ts, trans_df_chunk)
            cooling_accum = chunk_cf if cooling_accum is None else cooling_accum + chunk_cf
        if cooling_accum is None:
            cooling_func = np.zeros(len(Ts))
        else:
            cooling_func = cooling_accum
    else:
        trans_chunks = list(trans_reader)
        if len(trans_chunks) == 0:
            cooling_func = np.zeros(len(Ts))
        else:
            with _executor_context(max_workers=ncputrans) as trans_executor:
                futures = [trans_executor.submit(process_exomol_cooling_func_chunk, states_df, Ts, chunk)
                           for chunk in log_tqdm(trans_chunks, desc=desc)]
                cooling_func = np.sum([future.result() for future in log_tqdm(futures, desc='Combining '+trans_filename)], axis=0)
    return cooling_func

def save_exomol_cooling_func(states_df, Ntemp, Tmax):
    """
    Main function to calculate and save cooling functions for ExoMol database.

    Processes all transition files, calculates cooling functions at specified
    temperature intervals, and saves results in .cf format file.

    Parameters
    ----------
    states_df : pd.DataFrame
        Complete states DataFrame with 'id', 'E', 'g' columns
    Ntemp : int
        Temperature step interval
    Tmax : int
        Maximum temperature in Kelvin
    """
    # tqdm.write('Calculate cooling functions.') 
    print('Calculate cooling functions.')  
    # Import legacy-style configuration variables from core (set via Config.to_globals()).
    from pyexocross.core import read_path, data_info, save_path, ncpufiles, ncputrans  

    t = Timer()
    t.start()
    Ts = np.array(range(Ntemp, Tmax+1, Ntemp)) 
    print('Reading all transitions and calculating cooling functions ...')
    trans_filepaths = get_transfiles(read_path, data_info)
    # Process multiple files in parallel
    with _executor_context(max_workers=ncpufiles) as executor:
        # Submit reading tasks for each file
        futures = [executor.submit(process_exomol_cooling_func, states_df, Ts,
                                   trans_filepath) for trans_filepath in log_tqdm(trans_filepaths, desc='Cooling files')]
        cooling_func = sum(np.array([future.result() for future in futures]))
    
    cooling_func_df = pd.DataFrame()
    cooling_func_df['T'] = Ts
    cooling_func_df['cooling function'] = cooling_func
    t.end()
    print('Finished reading all transitions and calculating cooling functions!\n')
    
    print('Saving cooling functions into file ...')   
    ts = Timer()    
    ts.start()     
    cf_folder = save_path + 'cooling/'
    ensure_dir(cf_folder)
    cf_path = cf_folder + '__'.join(data_info[-2:]) + '.cf' 
    np.savetxt(cf_path, cooling_func_df, fmt="%8.1f %20.8E")
    ts.end()
    print_file_info('Cooling functions', ['T', 'Cooling function'], ['%8.1f','%20.8E'])
    print('Cooling functions file has been saved:', cf_path, '\n')  
    print('Cooling functions have been saved!\n')  
    # tqdm.write('Cooling functions has been saved!\n') 
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')  
