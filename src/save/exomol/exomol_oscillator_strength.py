"""
Save ExoMol oscillator strengths.

This module provides functions for calculating and saving oscillator strengths.
"""
import os
import numpy as np
import pandas as pd
import dask.dataframe as dd
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from src.base.utils import Timer, ensure_dir
from src.base.log import log_tqdm, print_file_info
from src.base.large_file import (
    is_large_trans_file,
    read_trans_chunks,
    save_large_txt,
)
from src.calculation.calculate_para import cal_v
from src.calculation.calculate_oscillator_strength import cal_oscillator_strength
from src.database.load_exomol import get_transfiles
from src.plot.plot_oscillator_strength import plot_oscillator_strength

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

# Process Oscillator Strength
def process_exomol_oscillator_strength_chunk(states_df,trans_df):
    """
    Process a chunk of transitions to calculate oscillator strengths.

    Merges transitions with state data and calculates oscillator strength (gf or f).

    Parameters
    ----------
    states_df : pd.DataFrame
        States DataFrame with 'id', 'E', 'g' columns
    trans_df : pd.DataFrame
        Transition DataFrame chunk with columns ['u', 'l', 'A']

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: u, l, os (oscillator strength), v (wavenumber)
    """
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
        return pd.DataFrame(columns=['u','l','os','v'])
    
    # Use vectorized lookup instead of merge (much faster and less memory)
    Ep = states_indexed.loc[trans_df['u'], 'E'].values
    Epp = states_indexed.loc[trans_df['l'], 'E'].values
    gp = states_indexed.loc[trans_df['u'], 'g'].values
    gpp = states_indexed.loc[trans_df['l'], 'g'].values
    A = trans_df['A'].values
    
    v = cal_v(Ep, Epp)
    # Filter out transitions where Ep < Epp (v < 0), skip invalid line list entries
    valid_v_mask = v >= 0
    if not valid_v_mask.any():
        return pd.DataFrame(columns=['u','l','os','v'])
    
    # Apply filter to all arrays
    gp = gp[valid_v_mask]
    gpp = gpp[valid_v_mask]
    A = A[valid_v_mask]
    v = v[valid_v_mask]
    u_values = trans_df['u'].values[valid_v_mask]
    l_values = trans_df['l'].values[valid_v_mask]
    
    oscillator_strength = cal_oscillator_strength(gp, gpp, A, v)
    
    oscillator_strength_df = pd.DataFrame({
        'u': u_values,
        'l': l_values,
        'os': oscillator_strength,
        'v': v
    })
    return (oscillator_strength_df)

def process_exomol_oscillator_strength(states_df, trans_filepath):
    """
    Process a single transition file to calculate oscillator strengths.

    Handles both large and small files, using appropriate processing strategies.

    Parameters
    ----------
    states_df : pd.DataFrame
        States DataFrame with 'id', 'E', 'g' columns
    trans_filepath : str
        Path to the transition file to process

    Returns
    -------
    pd.DataFrame
        Combined oscillator strength DataFrame from all chunks
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
        result_frames = []
        for trans_df_chunk in tqdm(trans_reader, desc=desc):
            chunk_df = process_exomol_oscillator_strength_chunk(states_df, trans_df_chunk)
            if not chunk_df.empty:
                result_frames.append(chunk_df)
        if result_frames:
            oscillator_strength_df = pd.concat(result_frames, ignore_index=True)
        else:
            oscillator_strength_df = pd.DataFrame(columns=['u','l','os','v'])
    else:
        trans_chunks = list(trans_reader)
        if len(trans_chunks) == 0:
            oscillator_strength_df = pd.DataFrame(columns=['u','l','os','v'])
        else:
            with _executor_context(max_workers=ncputrans) as trans_executor:
                futures = [trans_executor.submit(process_exomol_oscillator_strength_chunk, states_df, chunk)
                           for chunk in tqdm(trans_chunks, desc=desc)]
                oscillator_strength_df = pd.concat([future.result() for future in tqdm(futures, desc='Combining '+trans_filename)])
    return oscillator_strength_df

def save_exomol_oscillator_strength(states_df):
    """
    Main function to calculate and save oscillator strengths for ExoMol database.

    Processes all transition files, calculates oscillator strengths, and saves
    results in .os format file.

    Parameters
    ----------
    states_df : pd.DataFrame
        Complete states DataFrame with 'id', 'E', 'g' columns
    """
    # Import legacy-style configuration variables from core (set via Config.to_globals()).
    from pyexocross.core import ( 
        read_path,
        data_info,
        ncpufiles,
        save_path,
        database,
        gfORf,
        PlotOscillatorStrengthYN,
    )

    print('Calculate oscillator strengths.')  
    tot = Timer()
    tot.start()
    print('Reading transitions and calculating oscillator strengths ...')    
    trans_filepaths = get_transfiles(read_path, data_info)
    # Process multiple files in parallel
    with _executor_context(max_workers=ncpufiles) as executor:
        # Submit reading tasks for each file
        futures = [executor.submit(process_exomol_oscillator_strength,states_df,
                                   trans_filepath) for trans_filepath in tqdm(trans_filepaths)]
        oscillator_strength_df = pd.concat([future.result() for future in futures])
    oscillator_strength_df.sort_values(by=['v'], ascending=True, inplace=True)  
    # oscillator_strength_df  = oscillator_strength_df.sort_values('v').reset_index(drop=True)
    tot.end()
    print('Finished reading all transitions and calculating oscillator strengths!\n')

    print('Saving oscillator strengths into file ...')   
    ts = Timer()    
    ts.start()
    os_folder = save_path + 'oscillator_strength/files/'+data_info[0]+'/'+database+'/'
    ensure_dir(os_folder)
    os_path = os_folder+'__'.join(data_info[-2:])+'_'+gfORf.lower()+'.os' 
    os_format = "%12d %12d %10.4E %15.6f"
    save_large_txt(os_path, oscillator_strength_df, fmt=os_format)
    ts.end()
    print_file_info('Oscillator strengths', oscillator_strength_df.columns, ['%12d', '%12d', '%10.4E', '%15.6f'])
    print('Wavenumber in unit of cm⁻¹')
    print('Oscillator strengths file has been saved:', os_path, '\n')  
    
    # Plot oscillator strengths and save it as .png.
    if PlotOscillatorStrengthYN == 'Y':
        plot_oscillator_strength(oscillator_strength_df)
    print('\nOscillator strengths have been saved!\n')  
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')
