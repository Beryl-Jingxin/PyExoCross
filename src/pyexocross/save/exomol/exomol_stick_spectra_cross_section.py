# Combined function: Stick spectra and cross sections for ExoMol database
# Reads transition files ONCE for all temperatures and pressures
import numpy as np
import pandas as pd
from typing import TYPE_CHECKING
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

from pyexocross.base.config_manager import get_config
from pyexocross.base.utils import Timer, ensure_dir
from pyexocross.base.large_file import save_large_txt
from pyexocross.base.log import (
    log_tqdm, 
    print_stick_info, 
    print_xsec_info, 
    print_file_info,
    print_T_Tvib_Trot_P_path_info,
)
from pyexocross.base.large_file import (
    is_large_trans_file,
    read_trans_chunks,
    save_large_txt,
)
from pyexocross.database import read_broad, get_part_transfiles
from pyexocross.process.T_n_val import get_ntemp, get_temp_vals
from pyexocross.process.stick_xsec_filepath import stick_spectra_filepath, temperature_string_base
from pyexocross.save.exomol.exomol_stick_spectra import process_exomol_stick_spectra_chunk
from pyexocross.save.exomol.exomol_cross_section import process_exomol_cross_section_chunk
from pyexocross.plot.plot_stick_spectra import plot_stick_spectra
from pyexocross.plot.plot_cross_section import save_xsec_file_plot
from pyexocross.calculation.calcualte_line_profile import line_profile

if TYPE_CHECKING:
    from pyexocross.core import (  
        read_path,
        save_path,
        data_info,
        QNs_format,
        min_wnl,
        max_wnl,
        UncFilter,
        threshold,
        abs_emi,
        LTE_NLTE,
        photo,
        wn_wl,
        wn_wl_unit,
        PlotStickSpectraYN,
        PlotStickSpectraWnWl,
        PlotStickSpectraUnit,
        cutoff,
        wn_grid,
        alpha_hwhm_colid,
        gamma_hwhm_colid,
        NLTEMethod,
        ncputrans,
    )

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


def _ensure_exomol_globals():
    """
    Ensure legacy-style configuration globals are available in this module.
    """
    from pyexocross.core import ( 
        read_path,
        save_path,
        data_info,
        QNs_format,
        min_wnl,
        max_wnl,
        UncFilter,
        threshold,
        abs_emi,
        LTE_NLTE,
        photo,
        wn_wl,
        wn_wl_unit,
        PlotStickSpectraYN,
        PlotStickSpectraWnWl,
        PlotStickSpectraUnit,
        cutoff,
        wn_grid,
        alpha_hwhm_colid,
        gamma_hwhm_colid,
        NLTEMethod,
        ncputrans,
    )
    globals().update(
        dict(
            read_path=read_path,
            save_path=save_path,
            data_info=data_info,
            QNs_format=QNs_format,
            min_wnl=min_wnl,
            max_wnl=max_wnl,
            UncFilter=UncFilter,
            threshold=threshold,
            abs_emi=abs_emi,
            LTE_NLTE=LTE_NLTE,
            photo=photo,
            wn_wl=wn_wl,
            wn_wl_unit=wn_wl_unit,
            PlotStickSpectraYN=PlotStickSpectraYN,
            PlotStickSpectraWnWl=PlotStickSpectraWnWl,
            PlotStickSpectraUnit=PlotStickSpectraUnit,
            cutoff=cutoff,
            wn_grid=wn_grid,
            alpha_hwhm_colid=alpha_hwhm_colid,
            gamma_hwhm_colid=gamma_hwhm_colid,
            NLTEMethod=NLTEMethod,
            ncputrans=ncputrans,
        )
    )


def save_exomol_stick_spectra_cross_section(
    states_part_df,
    T_list,
    Tvib_list,
    Trot_list,
    P_list,
    Q_arr,
    check_uncertainty,
    check_lifetime,
    check_gfactor,
):
    """
    Combined function to calculate and save both stick spectra and cross sections.

    Reads transition files once and processes them for all temperature and pressure
    combinations to optimize memory usage.

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
    P_list : float, list, or list of lists
        Pressure(s) - single value, flat list, or list of lists per temperature
    Q_arr : np.ndarray
        Partition function array, shape (n_temps,)
    check_uncertainty : bool
        Whether states file contains uncertainty column and related handling is needed.
    check_lifetime : bool
        Whether states file contains lifetime column and related handling is needed.
    check_gfactor : bool
        Whether states file contains g-factor column and related handling is needed.
    """
    _ensure_exomol_globals()

    # def _format_temp_label(T_val, Tvib_val, Trot_val):
    #     temp_part = temperature_string_base(T_val, Tvib_val, Trot_val, NLTEMethod)
    #     if temp_part.startswith('Tvib'):
    #         temp_label = temp_part.replace('Tvib', 'Tvib=').replace('Trot', 'Trot=')
    #     else:
    #         temp_label = f'T={temp_part[1:]}' if temp_part.startswith('T') else temp_part
    #         temp_label = temp_label.replace('Trot', 'Trot=')
    #     return temp_label.replace('__', ' and ')
    # Load required settings from cached configuration (set by Config)
    (
        DopplerHWHMYN,
        LorentzianHWHMYN,
        database,
        predissocYN,
        check_predissoc,
        profile,
    ) = get_config()
    print('Calculate stick spectra and cross sections.\n')
    # Separate timers for stick spectra and cross sections.
    t_ss = Timer()
    t_xsec = Timer()
    
    # Prepare states data
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
    
    # Read broadening data for cross sections
    broad, ratio, nbroad, broad_dfs = read_broad(read_path)
    profile_label = line_profile(profile)
    
    print('Reading transitions ONCE for all temperatures (will be used for both stick spectra and cross sections) ...')    
    trans_filepaths = get_part_transfiles(read_path, data_info)
    
    # Prepare file format strings
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
    
    # Determine columns needed for cross sections
    if DopplerHWHMYN == 'U' and LorentzianHWHMYN == 'U':
        use_cols_xsec = [0,1,2,alpha_hwhm_colid, gamma_hwhm_colid]
        use_names_xsec = ['u','l','A','alpha_hwhm', 'gamma_hwhm']
    elif DopplerHWHMYN == 'U' and LorentzianHWHMYN != 'U':
        use_cols_xsec = [0,1,2,alpha_hwhm_colid]
        use_names_xsec = ['u','l','A','alpha_hwhm']
    elif DopplerHWHMYN != 'U' and LorentzianHWHMYN == 'U':
        use_cols_xsec = [0,1,2,gamma_hwhm_colid]
        use_names_xsec = ['u','l','A','gamma_hwhm']
    else:
        use_cols_xsec = [0,1,2]
        use_names_xsec = ['u','l','A']
    
    # Check if profile is pressure-dependent (Lorentzian/Voigt need pressure)
    pressure_dependent = profile_label not in ['Gaussian', 'Doppler']
    
    # ntemp: L=len(T_list), T/D=len(Trot_list), P=1 (no temperature dimension)
    ntemp = get_ntemp(NLTEMethod, T_list, Trot_list)
    # Determine pressure structure
    # P_list can be: single value, flat list (same pressures for all T), or list of lists (pressures per T)
    if isinstance(P_list, (list, np.ndarray)) and len(P_list) > 0:
        first_item = P_list[0]
        if isinstance(first_item, (list, np.ndarray)):
            # Multiple pressures per temperature: P_list is list of lists
            P_per_temp = [list(P_list[i]) if i < len(P_list) else list(P_list[0]) for i in range(ntemp)]
        else:
            # Flat list of pressures: apply the same pressure list to every temperature (cartesian product)
            P_values = list(P_list)
            P_per_temp = [P_values[:] for _ in range(ntemp)]
    else:
        # Single pressure for all temperatures
        P_single = P_list if isinstance(P_list, (int, float)) else (P_list[0] if len(P_list) > 0 else 1.0)
        P_per_temp = [[P_single] for _ in range(ntemp)]
    
    # Initialize results storage
    # Stick spectra: once per temperature (pressure doesn't matter)
    stick_spectra_results = {temp_idx: [] for temp_idx in range(ntemp)}
    # Cross sections: for pressure-dependent profiles, we'll process per (T,P) combination
    # For non-pressure-dependent, process once per temperature
    if pressure_dependent:
        # Store as dict: (temp_idx, press_idx) -> xsec array
        xsec_results = {}
    else:
        # Store as dict: temp_idx -> xsec array
        xsec_results = {temp_idx: None for temp_idx in range(ntemp)}
    
    # Cache transition chunks for small files (read once, use for all temperatures)
    trans_chunks_cache_ss = {}  # Cache for stick spectra chunks
    trans_chunks_cache_xsec = {}  # Cache for cross section chunks
    large_files_list = []
    
    # Determine which files are large and cache small file chunks
    for trans_filepath in trans_filepaths:
        trans_filename = trans_filepath.split('/')[-1]
        large_file = is_large_trans_file(trans_filepath)
        large_files_list.append((trans_filepath, large_file))
        
        if not large_file:
            # Cache chunks for small files
            use_cols_ss = [0,1,2]
            use_names_ss = ['u','l','A']
            trans_reader_ss = read_trans_chunks(trans_filepath, use_cols_ss, use_names_ss)
            trans_chunks_cache_ss[trans_filepath] = list(trans_reader_ss)
            num_chunks_ss = len(trans_chunks_cache_ss[trans_filepath])
            
            # Cache cross section chunks if different from stick spectra
            if use_cols_ss != use_cols_xsec:
                trans_reader_xsec = read_trans_chunks(trans_filepath, use_cols_xsec, use_names_xsec)
                trans_chunks_cache_xsec[trans_filepath] = list(trans_reader_xsec)
            else:
                trans_chunks_cache_xsec[trans_filepath] = trans_chunks_cache_ss[trans_filepath]
            num_chunks_xsec = len(trans_chunks_cache_xsec[trans_filepath])
            
            if num_chunks_ss == 0:
                print(f'Warning: No chunks found when caching small file {trans_filename}')
            else:
                print(f'Cached {num_chunks_ss} chunk(s) for stick spectra and {num_chunks_xsec} chunk(s) for cross sections from {trans_filename}')
    
    # Process all transition files ONCE
    print(f'Processing {len(large_files_list)} transition file(s) for {ntemp} temperature(s)...')
    if len(large_files_list) == 0:
        print('Warning: No transition files found to process.')
    else:
        for trans_filepath, large_file in log_tqdm(large_files_list, desc='Processing transition files'):
            trans_filename = trans_filepath.split('/')[-1]
            
            if large_file:
                # For large files: stream ONCE, process each chunk for ALL temperatures
                print(f'Streaming large file: {trans_filename} (processing for all temperatures)')
                trans_reader = read_trans_chunks(trans_filepath, use_cols_xsec, use_names_xsec)
                chunk_count = 0
                
                for trans_df_chunk in trans_reader:
                    chunk_count += 1
                    # Process this chunk for ALL temperatures
                    for temp_idx in range(ntemp):
                        T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
                        P = P_list[temp_idx] if isinstance(P_list, (list, np.ndarray)) and len(P_list) > temp_idx else P_list
                        
                        # Process stick spectra for this temperature
                        chunk_ss = process_exomol_stick_spectra_chunk(states_part_df_ss, T_list, Tvib_list, Trot_list, Q_arr, trans_df_chunk, temp_idx)
                        if chunk_ss is not None and len(chunk_ss) > 0:
                            stick_spectra_results[temp_idx].append(chunk_ss)
                        
                        # Process cross sections for this temperature
                        # For pressure-dependent profiles, we'll process separately per pressure later
                        # For non-pressure-dependent profiles, process now
                        if not pressure_dependent:
                            chunk_xsec = process_exomol_cross_section_chunk(states_part_df,T_list,Tvib_list,Trot_list,P,Q_arr,
                                                                   broad,ratio,nbroad,broad_dfs,profile_label,trans_df_chunk,temp_idx)
                            if chunk_xsec is not None:
                                if xsec_results[temp_idx] is None:
                                    xsec_results[temp_idx] = chunk_xsec.copy()
                                else:
                                    xsec_results[temp_idx] = xsec_results[temp_idx] + chunk_xsec
                
                if chunk_count == 0:
                    print(f'Warning: No chunks found in large file {trans_filename}')
            else:
                # For small files: use cached chunks, process for ALL temperatures
                trans_chunks_ss = trans_chunks_cache_ss.get(trans_filepath, [])
                trans_chunks_xsec = trans_chunks_cache_xsec.get(trans_filepath, [])
                
                if len(trans_chunks_ss) == 0 and len(trans_chunks_xsec) == 0:
                    print(f'Warning: No cached chunks found for file {trans_filename}')
                    continue
                
                # Process all temperatures for this file
                for temp_idx in range(ntemp):
                    T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
                    P = P_list[temp_idx] if isinstance(P_list, (list, np.ndarray)) and len(P_list) > temp_idx else P_list
                    
                    # Process stick spectra
                    if len(trans_chunks_ss) > 0:
                        with _executor_context(max_workers=ncputrans) as executor_ss:
                            futures_ss = [executor_ss.submit(process_exomol_stick_spectra_chunk, states_part_df_ss, T_list, Tvib_list, Trot_list, Q_arr, chunk, temp_idx)
                                         for chunk in trans_chunks_ss]
                            chunk_results_ss = [future.result() for future in futures_ss]
                            chunk_results_ss = [df for df in chunk_results_ss if df is not None and len(df) > 0]
                            if chunk_results_ss:
                                stick_spectra_results[temp_idx].extend(chunk_results_ss)
                    
                    # Process cross sections (only for non-pressure-dependent profiles)
                    if not pressure_dependent and len(trans_chunks_xsec) > 0:
                        with _executor_context(max_workers=ncputrans) as executor_xsec:
                            futures_xsec = [
                                executor_xsec.submit(process_exomol_cross_section_chunk,states_part_df,T_list,Tvib_list,Trot_list,P,Q_arr,
                                                     broad,ratio,nbroad,broad_dfs,profile_label,chunk,temp_idx) 
                                for chunk in trans_chunks_xsec
                            ]
                            chunk_results_xsec = [future.result() for future in futures_xsec]
                            for chunk_xsec in chunk_results_xsec:
                                if chunk_xsec is not None:
                                    if xsec_results[temp_idx] is None:
                                        xsec_results[temp_idx] = chunk_xsec.copy()
                                    else:
                                        xsec_results[temp_idx] = xsec_results[temp_idx] + chunk_xsec
    
    # Save stick spectra FIRST (before cross sections)
    # Combine and save results
    # For stick spectra: save num_T files (one per temperature, pressure doesn't matter)a
    t_ss.start()
    any_results_ss = False
    ss_file_count = 0
    print('\nCalculate stick spectra.')
    print_stick_info('cm⁻¹', 'cm/molecule')
    ss_columns = None
    ss_fmt_list = ss_fmt.split()
    for temp_idx in log_tqdm(range(ntemp), desc='\nProcessing stick spectra'):
        T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
        
        # Combine and save stick spectra results for this temperature (one file per temperature)
        if len(stick_spectra_results[temp_idx]) > 0:
            any_results_ss = True
            stick_spectra_df = pd.concat(stick_spectra_results[temp_idx], ignore_index=True)
            
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
            
            # Save stick spectra file for this temperature (shared naming)
            ss_path = stick_spectra_filepath(ss_folder, T, Tvib, Trot, str_min_wnl, str_max_wnl, unit_fn,
                                             data_info, wn_wl, UncFilter, threshold, database, abs_emi, LTE_NLTE, photo,
                                             NLTEMethod)
            
            ts = Timer()    
            ts.start()
            save_large_txt(ss_path, stick_spectra_df, fmt=ss_fmt)
            ss_file_count += 1
            if ss_columns is None:
                ss_columns = list(stick_spectra_df.columns)
            ts.end()
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
            stick_spectra_results[temp_idx] = []  # Clear the list
    
    if any_results_ss:
        if ss_columns is not None:
            print_file_info('Stick spectra', ss_columns, ss_fmt_list)
        # End and report stick spectra timing
        print('\nTotal running time for stick spectra:')
        t_ss.end()
        print(f'\nAll {ss_file_count} stick spectra files have been saved!\n')
        print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')
    
    # Print cross section info once at the beginning
    # Collect all unique pressures for display
    all_pressures = []
    for P_temp in P_per_temp:
        all_pressures.extend(P_temp)
    unique_pressures = sorted(list(set(all_pressures)))
    
    # Print cross section information once before processing
    print('Calculate cross sections.')
    print_xsec_info(profile_label, cutoff, UncFilter, min_wnl, max_wnl, 
                    'cm⁻¹', 'cm⁻¹/(molecule cm⁻²)', broad, ratio)
    
    # Process cross sections separately for each (T, P) combination to save memory
    # This is especially important for pressure-dependent profiles (Lorentzian/Voigt)
    # Loop structure: Outer loop over pressures, inner loop over temperatures
    # This gives: T1P1, T2P1, ..., TnP1, T1P2, T2P2, ..., TnP2, ..., T1Pn, T2Pn, ..., TnPn
    any_results_xsec = False
    xsec_file_count = 0
    # Start cross-section timer just before cross-section loops
    t_xsec.start()
    if pressure_dependent:
        # Get all unique pressures (sorted for consistent ordering)
        # Since P_per_temp might have different pressures per temperature, we need to collect all unique ones
        all_pressures_set = set()
        for P_temp in P_per_temp:
            all_pressures_set.update(P_temp)
        all_pressures_sorted = sorted(list(all_pressures_set))
        
        # Create list of all (T, P) combinations in the correct order: P first, then T
        tp_combinations = []
        for P in all_pressures_sorted:
            for temp_idx in range(ntemp):
                # Check if this pressure is valid for this temperature
                if P in P_per_temp[temp_idx]:
                    tp_combinations.append((temp_idx, P))
        
        total_combinations = len(tp_combinations)
        print(f'Processing cross sections for {ntemp} temperature(s) and {len(all_pressures_sorted)} pressure(s) ({total_combinations} total combinations)...')
        print(f'Order: For each pressure, process all temperatures (P1: T1-T{ntemp}, P2: T1-T{ntemp}, ...)')
        
        for combo_idx, (temp_idx, P) in enumerate(log_tqdm(tp_combinations, desc='\nProcessing cross sections')):
            T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
            
            # Process this (T, P) combination
            xsec = np.zeros_like(wn_grid, dtype=np.float64)  # Initialize with zeros instead of None
            
            # Process all transition files for this (T, P) combination
            for trans_filepath, large_file in large_files_list:
                trans_filename = trans_filepath.split('/')[-1]
                
                if large_file:
                    # Stream large files again for this (T, P) combination
                    trans_reader = read_trans_chunks(trans_filepath, use_cols_xsec, use_names_xsec)
                    for trans_df_chunk in trans_reader:
                        chunk_xsec = process_exomol_cross_section_chunk(states_part_df,T_list,Tvib_list,Trot_list,P,Q_arr,
                                                               broad,ratio,nbroad,broad_dfs,profile_label,trans_df_chunk,temp_idx)
                        if chunk_xsec is not None:
                            xsec = chunk_xsec if xsec is None else xsec + chunk_xsec
                else:
                    # Use cached chunks for small files
                    trans_chunks_xsec = trans_chunks_cache_xsec.get(trans_filepath, [])
                    if len(trans_chunks_xsec) > 0:
                        with _executor_context(max_workers=ncputrans) as executor_xsec:
                            futures_xsec = [
                                executor_xsec.submit(process_exomol_cross_section_chunk,states_part_df,T_list,Tvib_list,Trot_list,P,Q_arr,
                                                     broad,ratio,nbroad,broad_dfs,profile_label,chunk,temp_idx) 
                                for chunk in trans_chunks_xsec
                            ]
                            chunk_results_xsec = [future.result() for future in futures_xsec]
                            for chunk_xsec in chunk_results_xsec:
                                if chunk_xsec is not None:
                                    xsec = chunk_xsec if xsec is None else xsec + chunk_xsec
            
            # Save cross section for this (T, P) combination
            if xsec is not None and len(xsec) > 0 and not np.all(xsec == 0):
                any_results_xsec = True
                # Create a unique identifier for this (T, P) combination
                print_T_Tvib_Trot_P_path_info(T, Tvib, Trot, P, abs_emi, NLTEMethod, 'Cross sections', None)
                save_xsec_file_plot(wn_grid, xsec, database, profile_label, T, P, temp_idx, Tvib_list, Trot_list)
                xsec_file_count += 1
                del xsec
    
    # Process cross section results (only for non-pressure-dependent profiles)
    # Pressure-dependent profiles are already processed and saved above
    if not pressure_dependent:
        any_results_xsec = False
        for temp_idx in log_tqdm(range(ntemp), desc='\nSaving cross sections'):
            T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
            if xsec_results.get(temp_idx) is not None and len(xsec_results[temp_idx]) > 0 and not np.all(xsec_results[temp_idx] == 0):
                # For non-pressure-dependent profiles, save one file per temperature
                any_results_xsec = True
                P = P_per_temp[temp_idx][0]  # Use first pressure for non-pressure-dependent
                # Save cross section for this temperature
                print_T_Tvib_Trot_P_path_info(T, Tvib, Trot, None, abs_emi, NLTEMethod, 'Cross sections', None)
                save_xsec_file_plot(wn_grid, xsec_results[temp_idx], database, profile_label, T, P, temp_idx, Tvib_list, Trot_list)
                xsec_file_count += 1
                # Clear memory
                xsec_results[temp_idx] = None
    
    if any_results_xsec:
        print_file_info('Cross sections', ['Wavenumber', 'Cross section'], ['%12.6f', '%14.8E'])
        print('\nTotal running time for cross sections:')
        t_xsec.end()
        print(f'\nAll {xsec_file_count} cross sections files have been saved!\n')
    print('\nFinished reading transitions and calculating stick spectra and cross sections!\n')
    
    # Check if we have any results
    total_ss_results = ss_file_count
    if pressure_dependent:
        total_xsec_results = xsec_file_count
    else:
        total_xsec_results = xsec_file_count
    
    if not any_results_ss and not any_results_xsec:
        print(f'Info: Total stick spectra results accumulated: {total_ss_results}')
        print(f'Info: Total cross section results accumulated: {total_xsec_results}')
        print(f'Info: Number of transition files processed: {len(large_files_list)}')
        raise ValueError("Empty result with the input filter values. Please type new filter values in the input file.")
    
    if any_results_ss:
        print(f'All {ss_file_count} stick spectra files have been saved!\n')
    if any_results_xsec:
        print(f'All {xsec_file_count} cross sections files have been saved!\n')
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')
