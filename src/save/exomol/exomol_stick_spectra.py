"""
Save ExoMol stick spectra.

This module provides functions for calculating and saving stick spectra
from ExoMol database files.
"""
import os
from functools import partial
import numpy as np
import pandas as pd
import dask.dataframe as dd
from concurrent.futures import ProcessPoolExecutor
from src.base.utils import Timer, ensure_dir
from src.base.log import log_tqdm, print_stick_info, print_T_Tvib_Trot_P_path_info
from src.base.large_file import (
    save_large_txt,
    is_large_trans_file,
    read_trans_chunks,
    process_large_chunks,
)
from src.base.constants import MAX_LARGE_FILE_WORKERS
from src.calculation.calculate_para import cal_v
from src.calculation.calculate_intensity import (
    cal_abscoefs,
    cal_abscoefs_nlte_2T,
    cal_abscoefs_nlte_nvib,
    cal_abscoefs_nlte_pop,
)
from src.calculation.calculate_emissivity import (
    cal_emicoefs,
    cal_emicoefs_nlte_2T,
    cal_emicoefs_nlte_nvib,
    cal_emicoefs_nlte_pop,
)
from src.process.T_n_val import get_ntemp, get_temp_vals
from src.process.S_for_LTE_NLTE_Ab_Em import S_for_LTE_NLTE_Ab_Em
from src.process.stick_xsec_filepath import stick_spectra_filepath
from src.database import get_part_transfiles
from src.process.filter_qn import QNfilter_linelist
from src.plot.plot_stick_spectra import plot_stick_spectra

# Stick Spectra
# Process LTE or NLTE linelist
def process_exomol_stick_spectra_chunk(states_part_df,T_list,Tvib_list,Trot_list,Q_arr,trans_part_df,temp_idx=None):
    """
    Process a chunk of transitions to calculate stick spectra intensities.

    Merges transitions with state data, calculates wavenumbers, applies filters,
    and computes absorption or emission coefficients based on LTE or non-LTE method.

    Parameters
    ----------
    states_part_df : pd.DataFrame
        Filtered states DataFrame
    T_list : list of float
        Temperature list for LTE calculations
    Tvib_list : list of float
        Vibrational temperature list for non-LTE two-temperature method
    Trot_list : list of float
        Rotational temperature list for non-LTE methods
    Q_arr : np.ndarray
        Partition function array, shape (n_temps,)
    trans_part_df : pd.DataFrame
        Transition DataFrame chunk with columns ['u', 'l', 'A']
    temp_idx : int, optional
        Temperature index to process single temperature (for memory efficiency)

    Returns
    -------
    pd.DataFrame
        Stick spectra DataFrame with columns: v, S, J', E', J", E", and quantum numbers
    """
    from pyexocross.core import ( 
        min_wn,
        max_wn,
        QNsFilter,
        QNs_value,
        QNs_label,
        abs_emi,
        NLTEMethod,
        abundance,
        threshold,
    )
    # Optimized: use indexed lookup instead of two merges to reduce memory
    # If temp_idx is provided, process only that temperature to save memory
    if isinstance(states_part_df, dd.DataFrame):
        states_part_df = states_part_df.compute()
    if isinstance(trans_part_df, dd.DataFrame):
        trans_part_df = trans_part_df.compute()
    
    # Set index for fast lookup (more memory efficient than merge)
    states_indexed = states_part_df.set_index('id')
    
    # Filter trans_df to only include transitions where both states exist
    valid_mask = trans_part_df['u'].isin(states_indexed.index) & trans_part_df['l'].isin(states_indexed.index)
    trans_part_df = trans_part_df[valid_mask]
    
    if len(trans_part_df) == 0:
        states_cols = [col for col in states_part_df.columns if col != 'id']
        u_cols = [col + "'" for col in states_cols]
        l_cols = [col + '"' for col in states_cols]
        combined_cols = ['A'] + u_cols + l_cols
        col_main = ['v','S',"J'","E'",'J"','E"']
        col_qn = [col for col in combined_cols if col not in col_main]
        col_stick_spectra = col_main + col_qn
        return pd.DataFrame(columns=col_stick_spectra)
    
    # Get upper and lower state data using vectorized lookup
    u_states = states_indexed.loc[trans_part_df['u']].reset_index(drop=True)
    l_states = states_indexed.loc[trans_part_df['l']].reset_index(drop=True)
    
    # Rename columns with suffixes (id column is already dropped by reset_index)
    u_states.columns = [col + "'" for col in u_states.columns]
    l_states.columns = [col + '"' for col in l_states.columns]
    
    # Combine trans and states data
    stick_spectra_df = pd.concat([
        trans_part_df[['A']].reset_index(drop=True),
        u_states,
        l_states
    ], axis=1)
    
    stick_spectra_df['v'] = cal_v(stick_spectra_df["E'"].values, stick_spectra_df['E"'].values)
    # Filter out transitions where Ep < Epp (v < 0), skip invalid line list entries
    stick_spectra_df = stick_spectra_df[stick_spectra_df['v'] >= 0]
    stick_spectra_df = stick_spectra_df[stick_spectra_df['v'].between(min_wn, max_wn)]
    if len(stick_spectra_df) != 0 and QNsFilter != []:
        stick_spectra_df = QNfilter_linelist(stick_spectra_df, QNs_value, QNs_label)
    num = len(stick_spectra_df)
    if num > 0:
        stick_spectra_df['S'] = S_for_LTE_NLTE_Ab_Em(stick_spectra_df,T_list,Tvib_list,Trot_list,Q_arr,abs_emi,NLTEMethod,temp_idx)
        if threshold != 'None':
            stick_spectra_df = stick_spectra_df[stick_spectra_df['S'] >= threshold]
        else:
            pass
    else:
        stick_spectra_df['S'] = np.array([])

    if NLTEMethod == 'L':
        col_list = ['A',"g'",'g"']
    elif NLTEMethod == 'T':
        col_list = ['A',"g'",'g"',"Evib'","Erot'",'Evib"','Erot"']
    elif NLTEMethod == 'D':
        col_list = ['A',"g'",'g"',"nvib'",'nvib"']
    elif NLTEMethod == 'P':
        col_list = ['A',"g'",'g"',"pop'",'pop"']
    else:
        raise ValueError("Please choose one non-LTE method from: 'T', 'D' or 'P'.")
    stick_spectra_df.drop(columns=col_list, inplace=True)
    col_main = ['v','S',"J'","E'",'J"','E"']
    col_qn = [col for col in stick_spectra_df.columns if col not in col_main]
    col_stick_spectra = col_main + col_qn
    stick_spectra_df = stick_spectra_df[col_stick_spectra]  
    return stick_spectra_df

def process_exomol_stick_spectra(states_part_df,T_list,Tvib_list,Trot_list,Q_arr,trans_filepath,temp_idx=None):
    """
    Process a single transition file to calculate stick spectra.

    Handles both large and small files, using appropriate processing strategies
    for memory efficiency.

    Parameters
    ----------
    states_part_df : pd.DataFrame
        Filtered states DataFrame
    T_list : list of float
        Temperature list for LTE calculations
    Tvib_list : list of float
        Vibrational temperature list for non-LTE two-temperature method
    Trot_list : list of float
        Rotational temperature list for non-LTE methods
    Q_arr : np.ndarray
        Partition function array, shape (n_temps,)
    trans_filepath : str
        Path to the transition file to process
    temp_idx : int, optional
        Temperature index to process single temperature

    Returns
    -------
    pd.DataFrame
        Combined stick spectra DataFrame from all chunks
    """   
    from pyexocross.core import ncputrans
    trans_filename = trans_filepath.split('/')[-1]
    print('Processeing transitions file:', trans_filename)
    use_cols = [0,1,2]
    use_names = ['u','l','A']
    large_file = is_large_trans_file(trans_filepath)
    trans_reader = read_trans_chunks(trans_filepath, use_cols, use_names)
    desc = 'Processing ' + trans_filename + (' (limited streaming)' if large_file else '')
    zero_factory = lambda: pd.DataFrame(columns=['v','S',"J'","E'",'J"','E"'])
    def combine_fn(results):
        if not results:
            return zero_factory()
        return pd.concat(results, ignore_index=True)
    handler = partial(process_exomol_stick_spectra_chunk, states_part_df, T_list, Tvib_list, Trot_list, Q_arr, temp_idx=temp_idx)
    if large_file:
        print('Large transition file detected (>1 GB). Using bounded parallel streaming to reduce memory usage.')
        stick_spectra_df = process_large_chunks(
            trans_reader,
            handler,
            combine_fn,
            zero_factory,
            desc,
            max_workers=max(1, min(ncputrans, MAX_LARGE_FILE_WORKERS))
        )
    else:
        trans_chunks = list(trans_reader)
        if len(trans_chunks) == 0:
            stick_spectra_df = zero_factory()
        else:
            with ProcessPoolExecutor(max_workers=ncputrans) as trans_executor:
                futures = [trans_executor.submit(process_exomol_stick_spectra_chunk, states_part_df, T_list, Tvib_list, Trot_list, Q_arr, chunk, temp_idx)
                           for chunk in log_tqdm(trans_chunks, desc=desc)]
                stick_spectra_df = combine_fn([future.result() for future in log_tqdm(futures, desc='Combining '+trans_filename)])
    return stick_spectra_df

# Stick spectra for ExoMol database
def save_exomol_stick_spectra(states_part_df, T_list, Tvib_list, Trot_list, Q_arr):
    """
    Main function to calculate and save stick spectra for ExoMol database.

    Processes all transition files, calculates stick spectra for each temperature,
    and saves results in .stick format files.

    Parameters
    ----------
    states_part_df : pd.DataFrame
        Filtered states DataFrame
    T_list : list of float
        Temperature list for LTE calculations
    Tvib_list : list of float
        Vibrational temperature list for non-LTE two-temperature method
    Trot_list : list of float
        Rotational temperature list for non-LTE methods
    Q_arr : np.ndarray
        Partition function array, shape (n_temps,)
    """
    from pyexocross.core import (
        check_uncertainty,
        check_lifetime,
        check_gfactor,
        predissocYN,
        read_path,
        data_info,
        QNs_format,
        save_path,
        database,
        min_wnl,
        max_wnl,
        UncFilter,
        threshold,
        NLTEMethod,
        LTE_NLTE,
        photo,
        abs_emi,
        wn_wl,
        wn_wl_unit,
        PlotStickSpectraYN,
        ncpufiles,
        ncputrans,
    )
    print('Calculate stick spectra.')  
    print_stick_info('cm⁻¹', 'cm/molecule')
    tot = Timer()
    tot.start()
    # Q = read_exomol_pf(read_path, T)
    states_part_df_ss = states_part_df.copy()
    if check_uncertainty == True:
        states_part_df_ss.drop(columns=['unc'], inplace=True)
    if check_lifetime == True or predissocYN == 'Y':
        try:
            states_part_df_ss.drop(columns=['tau'], inplace=True)
        except:
            pass
    if check_gfactor == True:
        try:
            states_part_df_ss.drop(columns=['gfac'], inplace=True)
        except:
            pass
    
    print('\nReading transitions and calculating stick spectra ...')    
    trans_filepaths = get_part_transfiles(read_path, data_info)
    
    # Process each temperature separately to save memory
    QNsfmf = (str(QNs_format).replace("'","").replace(",","").replace("[","").replace("]","")
                .replace('d','s').replace('i','s').replace('.1f','s'))
    ss_folder = save_path + 'stick_spectra/files/'+data_info[0]+'/'+database+'/'
    ensure_dir(ss_folder)
    str_min_wnl = str(int(np.floor(min_wnl)))
    str_max_wnl = str(int(np.ceil(max_wnl)))
    
    if QNsfmf == '':
        ss_fmt = '%12.8E %12.8E %7s %12.4f %7s %12.4f'
    else:
        ss_fmt = '%12.8E %12.8E %7s %12.4f %7s %12.4f ' + QNsfmf + ' ' + QNsfmf
    
    # Process each temperature separately (ntemp: L=len(T_list), T/D=len(Trot_list), P=1)
    ntemp = get_ntemp(NLTEMethod, T_list, Trot_list)
    any_results = False
    ss_file_count = 0
    for temp_idx in log_tqdm(range(ntemp), desc='\nProcessing stick spectra'):
        T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
        
        # Process multiple files in parallel for this temperature
        with ProcessPoolExecutor(max_workers=ncpufiles) as executor:
            # Submit reading tasks for each file
            futures = [executor.submit(process_exomol_stick_spectra,states_part_df_ss,T_list,Tvib_list,Trot_list,Q_arr,
                                       trans_filepath,temp_idx) for trans_filepath in trans_filepaths]
            stick_spectra_df = pd.concat([future.result() for future in futures])
            
        if len(stick_spectra_df) == 0:
            print(f'Warning: No transitions found for T={T} K. Skipping this temperature.')
            continue
        
        any_results = True
            
        if wn_wl == 'WN':
            print_unit_str = 'Wavenumber in unit of cm⁻¹'
            unit_fn = 'cm-1__'
        elif wn_wl == 'WL' and wn_wl_unit == 'um':
            stick_spectra_df['v'] = 1e4/stick_spectra_df['v']
            print_unit_str = 'Wavelength in unit of μm'
            unit_fn = 'um__'
        elif wn_wl == 'WL' and wn_wl_unit == 'nm':
            stick_spectra_df['v'] = 1e7/stick_spectra_df['v']
            print_unit_str = 'Wavelength in unit of nm'
            unit_fn = 'nm__'
        else:
            raise ValueError('Please wirte the unit of wavelength in the input file: um or nm.')      
        stick_spectra_df.sort_values(by=['v'], ascending=True, inplace=True) 
        
        # Save file for this temperature (shared naming)
        ss_path = stick_spectra_filepath(ss_folder, T, Tvib, Trot, str_min_wnl, str_max_wnl, unit_fn,
                                        data_info, wn_wl, UncFilter, threshold, database, abs_emi, LTE_NLTE, photo,
                                        NLTEMethod)
        
        ts = Timer()    
        ts.start()
        save_large_txt(ss_path, stick_spectra_df, fmt=ss_fmt)
        ts.end()
        ss_file_count += 1
        print_T_Tvib_Trot_P_path_info(T, Tvib, Trot, None, abs_emi, NLTEMethod, 'Stick spectra', ss_path)

        # Plot stick spectra for this temperature
        if PlotStickSpectraYN == 'Y':
            # Create a copy for plotting
            stick_spectra_df_plot = stick_spectra_df.copy()
            # plot_stick_spectra expects wavenumber (cm⁻¹) input, so convert if needed
            if wn_wl == 'WL':
                # Convert back to wavenumber for plotting
                if wn_wl_unit == 'um':
                    stick_spectra_df_plot['v'] = 1e4/stick_spectra_df_plot['v']
                elif wn_wl_unit == 'nm':
                    stick_spectra_df_plot['v'] = 1e7/stick_spectra_df_plot['v']
            # If wn_wl == 'WN', data is already in wavenumber, so use it directly
            plot_stick_spectra(stick_spectra_df_plot, T=T, Tvib=Tvib, Trot=Trot)
            del stick_spectra_df_plot
        
        # Clear memory
        del stick_spectra_df
        
    tot.end()
    print('\nFinished reading transitions and calculating stick spectra!\n')
    
    if not any_results:
        raise ValueError("Empty result with the input filter values. Please type new filter values in the input file.")
    
    print(f'All {ss_file_count} stick spectra files have been saved!\n')
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')
