"""
Save ExoMol cross sections.

This module provides functions for calculating and saving cross sections
from ExoMol database files.
"""
import os
import numpy as np
import pandas as pd
import dask.dataframe as dd
from tqdm import tqdm
from functools import partial
from concurrent.futures import ProcessPoolExecutor
from src.base.utils import Timer
from src.base.log import log_tqdm, print_xsec_info, print_T_Tvib_Trot_P_path_info
from src.base.large_file import (
    is_large_trans_file,
    read_trans_chunks
)
from src.database import read_broad, get_part_transfiles
from src.database.load_exomol import extract_broad
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
from src.calculation.calcualte_line_profile import (
    DopplerHWHM_alpha,
    LorentzianHWHM_gamma,
    Gaussian_standard_deviation,
    line_profile,
    PseudoThompsonVoigt,
    PseudoKielkopfVoigt,
    PseudoOliveroVoigt,
    PseudoLiuLinVoigt,
    PseudoRoccoVoigt,
)
from src.process.filter_qn import QNfilter_linelist
from src.calculation.calculate_cross_section import (
    cross_section_Doppler,
    cross_section_Lorentzian,
    cross_section_SciPyVoigt,
    cross_section_SciPyWofzVoigt,
    cross_section_HumlicekVoigt,
    cross_section_PseudoVoigt,
    cross_section_BinnedGaussian,
    cross_section_BinnedLorentzian,
    cross_section_BinnedVoigt
)
from src.process.T_n_val import get_ntemp, get_temp_vals
from src.process.S_for_LTE_NLTE_Ab_Em import S_for_LTE_NLTE_Ab_Em
from src.plot.plot_cross_section import save_xsec_file_plot

## Line list for calculating cross sections
def process_exomol_cross_section_chunk(states_part_df,T_list,Tvib_list,Trot_list,P,Q_arr,
                                       broad,ratio,nbroad,broad_dfs,profile_label,
                                       trans_part_df,temp_idx=None):
    """
    Process a chunk of transitions to calculate cross sections.

    Merges transitions with state data, calculates line strengths, applies line profiles,
    and computes cross sections on the wavenumber grid.

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
    P : float
        Pressure value in bar
    Q_arr : np.ndarray
        Partition function array, shape (n_temps,)
    broad : list of str
        List of broadener names
    ratio : list of float
        Mixing ratios for broadeners
    nbroad : int
        Number of broadeners
    broad_dfs : list of pd.DataFrame
        Broadening DataFrames for each broadener
    profile_label : str
        Line profile label (e.g., 'Voigt', 'Doppler', 'Lorentzian')
    trans_part_df : pd.DataFrame
        Transition DataFrame chunk
    temp_idx : int, optional
        Temperature index to process single temperature

    Returns
    -------
    np.ndarray
        Cross-section array, shape (n_points,)
    """
    from pyexocross.core import (  
        wn_grid,
        cutoff,
        min_wn,
        max_wn,
        QNsFilter,
        QNs_value,
        QNs_label,
        NLTEMethod,
        abs_emi,
        abundance,
        predissocYN,
        profile,
        DopplerHWHMYN,
        LorentzianHWHMYN,
        threshold,
        alpha_hwhm_colid,
        gamma_hwhm_colid,
        alpha_HWHM as config_alpha_HWHM,
        gamma_HWHM as config_gamma_HWHM,
    )
    globals().update(
        dict(
            alpha_HWHM=config_alpha_HWHM,
            gamma_HWHM=config_gamma_HWHM,
        )
    )
    # Optimized: use indexed lookup instead of two merges to reduce memory
    # Note: P is a single pressure value, not a list
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
        return np.zeros_like(wn_grid)
    
    # Get upper and lower state data using vectorized lookup
    u_states = states_indexed.loc[trans_part_df['u']].reset_index(drop=True)
    l_states = states_indexed.loc[trans_part_df['l']].reset_index(drop=True)
    
    # Rename columns with suffixes (id column is already dropped by reset_index)
    u_states.columns = [col + "'" for col in u_states.columns]
    l_states.columns = [col + '"' for col in l_states.columns]
    
    # Combine trans and states data
    st_df = pd.concat([
        trans_part_df[['A']].reset_index(drop=True),
        u_states,
        l_states
    ], axis=1)
    
    st_df['v'] = cal_v(st_df["E'"].values, st_df['E"'].values)
    # Filter out transitions where Ep < Epp (v < 0), skip invalid line list entries
    st_df = st_df[st_df['v'] >= 0]
    if cutoff == 'None':
        st_df = st_df[st_df['v'].between(min_wn, max_wn)]
    else:
        st_df = st_df[st_df['v'].between(min_wn - cutoff, max_wn + cutoff)]
    if len(st_df) != 0 and QNsFilter != []:
        st_df = QNfilter_linelist(st_df, QNs_value, QNs_label)

    if NLTEMethod == 'L':
        extra_col = []
    elif NLTEMethod == 'T':
        if abs_emi == 'Ab':
            extra_col = ['Evib"', 'Erot"']
        else:
            extra_col = ["Evib'", "Erot'"]
    elif NLTEMethod == 'D':
        if abs_emi == 'Ab':
            extra_col = ['nvib"']
        else:
            extra_col = ["nvib'"]
    elif NLTEMethod == 'P':
        if abs_emi == 'Ab':
            extra_col = ['pop"', 'g"']
        else:
            extra_col = ["pop'"]
    else:
        extra_col = []

    if predissocYN == 'Y' and 'VOI' in profile:
        if DopplerHWHMYN == 'U' and LorentzianHWHMYN == 'U':
            st_df = st_df[['A','v',"g'","E'",'E"','J"',"tau'",'alpha_hwhm','gamma_hwhm']+extra_col]
        elif DopplerHWHMYN == 'U' and LorentzianHWHMYN != 'U':
            st_df = st_df[['A','v',"g'","E'",'E"','J"',"tau'",'alpha_hwhm']+extra_col]
        elif DopplerHWHMYN != 'U' and LorentzianHWHMYN == 'U':
            st_df = st_df[['A','v',"g'","E'",'E"','J"',"tau'",'gamma_hwhm']+extra_col]
        else:
            st_df = st_df[['A','v',"g'","E'",'E"','J"',"tau'"]+extra_col]
    else:
        if DopplerHWHMYN == 'U' and LorentzianHWHMYN == 'U':
            st_df = st_df[['A','v',"g'","E'",'E"','J"','alpha_hwhm','gamma_hwhm']+extra_col]
        elif DopplerHWHMYN == 'U' and LorentzianHWHMYN != 'U':
            st_df = st_df[['A','v',"g'","E'",'E"','J"','alpha_hwhm']+extra_col]
        elif DopplerHWHMYN != 'U' and LorentzianHWHMYN == 'U':
            st_df = st_df[['A','v',"g'","E'",'E"','J"','gamma_hwhm']+extra_col]
        else:
            st_df = st_df[['A','v',"g'","E'",'E"','J"']+extra_col]
        
    num = len(st_df)
    if num > 0:
        P_val = P
        T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
        st_df['coef'] = S_for_LTE_NLTE_Ab_Em(st_df,T_list,Tvib_list,Trot_list,Q_arr,abs_emi,NLTEMethod,temp_idx)
        st_df.drop(columns=['A',"g'","E'",'E"'], inplace=True)
        if threshold != 'None':
            st_df = st_df[st_df['coef'] >= threshold]  
        v = st_df['v'].values
        num = len(st_df)         

        global alpha_HWHM, gamma_HWHM
        if predissocYN == 'Y' and 'VOI' in profile:
            if DopplerHWHMYN == 'U' and LorentzianHWHMYN == 'U':
                alpha_HWHM = st_df['alpha_hwhm'].values
                gamma_HWHM = st_df['gamma_hwhm'].values
            elif DopplerHWHMYN == 'U' and LorentzianHWHMYN != 'U':
                alpha_HWHM = st_df['alpha_hwhm'].values
            elif DopplerHWHMYN != 'U' and LorentzianHWHMYN == 'U':
                gamma_HWHM = st_df['gamma_hwhm'].values
            else:
                pass
        else:
            if DopplerHWHMYN == 'U' and LorentzianHWHMYN != 'U':
                alpha_HWHM = st_df['alpha_hwhm'].values
            elif DopplerHWHMYN != 'U' and LorentzianHWHMYN == 'U':
                gamma_HWHM = st_df['gamma_hwhm'].values
            else:
                pass

        if num > 0:
            T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
            if NLTEMethod == 'T':
                if len(Trot_list) > 0:
                    trot_idx = temp_idx if temp_idx is not None and temp_idx < len(Trot_list) else 0
                    T = Trot_list[trot_idx]
            elif temp_idx is not None and temp_idx >= len(T_list):
                T = T_list[0] if len(T_list) > 0 else None
            if 'DOP' not in profile and 'GAU' not in profile:
                # Use dictionaries instead of DataFrames to avoid indexing issues
                gamma_L = {}
                n_air = {}
                for i in range(nbroad):
                    if broad[i] == 'Default':
                        gamma_default = float(broad_dfs[i]['gamma_L'][0])
                        n_air_default = float(broad_dfs[i]['n_air'][0])
                        gamma_L[i] = np.full(num, gamma_default, dtype=np.float64) * ratio[i]
                        n_air[i] = np.full(num, n_air_default, dtype=np.float64) * ratio[i]
                    else:
                        (gammaL,nair) = extract_broad(broad_dfs[i], st_df)
                        # extract_broad returns pandas Series, convert to numpy array
                        gamma_L[i] = np.asarray(gammaL, dtype=np.float64) * ratio[i]
                        n_air[i] = np.asarray(nair, dtype=np.float64) * ratio[i]
                        # Ensure correct shape
                        if len(gamma_L[i]) != num:
                            raise ValueError(f"gamma_L[{i}] from extract_broad has length {len(gamma_L[i])}, expected {num}")
                        if len(n_air[i]) != num:
                            raise ValueError(f"n_air[{i}] from extract_broad has length {len(n_air[i])}, expected {num}")    
                if predissocYN == 'Y' and 'VOI' in profile:
                    tau = st_df["tau'"].values
                else:
                    tau = np.array([])
                gamma = LorentzianHWHM_gamma(num, gamma_HWHM, nbroad, gamma_L, n_air, [], [], tau, T, P_val)
            else:
                gamma = np.array([])
            coef = st_df['coef'].values
        
            # Line profiles
            if profile_label == 'Doppler':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                xsec = cross_section_Doppler(wn_grid, v, alpha, coef, cutoff)
            elif profile_label == 'Gaussian':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                xsec = cross_section_Doppler(wn_grid, v, alpha, coef, cutoff)
            elif profile_label == 'Lorentzian':
                xsec = cross_section_Lorentzian(wn_grid, v, gamma, coef, cutoff)
            elif profile_label == 'SciPy Voigt':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                sigma = Gaussian_standard_deviation(alpha)     
                xsec = cross_section_SciPyVoigt(wn_grid, v, sigma, gamma, coef, cutoff)
            elif profile_label == 'SciPy wofz Voigt':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                sigma = Gaussian_standard_deviation(alpha)     
                xsec = cross_section_SciPyWofzVoigt(wn_grid, v, sigma, gamma, coef, cutoff)
            elif profile_label == 'Humlicek Voigt':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                xsec = cross_section_HumlicekVoigt(wn_grid, v, alpha, gamma, coef, cutoff)  
            elif profile_label == 'Thompson pseudo-Voigt':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                eta = PseudoThompsonVoigt(alpha, gamma)
                xsec = cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, coef, cutoff)       
            elif profile_label == 'Kielkopf pseudo-Voigt':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                eta = PseudoKielkopfVoigt(alpha, gamma)
                xsec = cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, coef, cutoff)       
            elif profile_label == 'Olivero pseudo-Voigt':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                eta = PseudoOliveroVoigt(alpha, gamma)
                xsec = cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, coef, cutoff)       
            elif profile_label == 'Liu-Lin pseudo-Voigt':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                eta = PseudoLiuLinVoigt(alpha, gamma)
                xsec = cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, coef, cutoff)       
            elif profile_label == 'Rocco pseudo-Voigt':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                eta = PseudoRoccoVoigt(alpha, gamma)
                xsec = cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, coef, cutoff) 
            elif profile_label == 'Binned Doppler':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                xsec = cross_section_BinnedGaussian(wn_grid, v, alpha, coef, cutoff)        
            elif profile_label == 'Binned Gaussion':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                xsec = cross_section_BinnedGaussian(wn_grid, v, alpha, coef, cutoff)
            elif profile_label == 'Binned Lorentzian':
                xsec = cross_section_BinnedLorentzian(wn_grid, v, gamma, coef, cutoff)
            elif profile_label == 'Binned Voigt':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                sigma = Gaussian_standard_deviation(alpha)     
                xsec = cross_section_BinnedVoigt(wn_grid, v, sigma, gamma, coef, cutoff)           
            else:
                raise ValueError('Please choose line profile from the list.')  
        else:      
                xsec = np.zeros_like(wn_grid)
    else:
        xsec = np.zeros_like(wn_grid)
    return xsec

def process_exomol_cross_section(states_part_df,T_list,Tvib_list,Trot_list,P,Q_arr,broad,ratio,nbroad,broad_dfs,profile_label,trans_filepath,temp_idx=None):
    """
    Process a single transition file to calculate cross sections.

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
    P : float
        Pressure value in bar
    Q_arr : np.ndarray
        Partition function array, shape (n_temps,)
    broad : list of str
        List of broadener names
    ratio : list of float
        Mixing ratios for broadeners
    nbroad : int
        Number of broadeners
    broad_dfs : list of pd.DataFrame
        Broadening DataFrames for each broadener
    profile_label : str
        Line profile label
    trans_filepath : str
        Path to the transition file to process
    temp_idx : int, optional
        Temperature index to process single temperature

    Returns
    -------
    np.ndarray
        Combined cross-section array from all chunks, shape (n_points,)
    """ 
    from pyexocross.core import (
        DopplerHWHMYN,
        LorentzianHWHMYN,
        alpha_hwhm_colid,
        gamma_hwhm_colid,
        wn_grid,
        ncputrans,
    )
    trans_filename = trans_filepath.split('/')[-1]
    print('Processeing transitions file:', trans_filename)
    if DopplerHWHMYN == 'U' and LorentzianHWHMYN == 'U':
        use_cols = [0,1,2,alpha_hwhm_colid, gamma_hwhm_colid]
        use_names = ['u','l','A','alpha_hwhm', 'gamma_hwhm']
    elif DopplerHWHMYN == 'U' and LorentzianHWHMYN != 'U':
        use_cols = [0,1,2,alpha_hwhm_colid]
        use_names = ['u','l','A','alpha_hwhm']
    elif DopplerHWHMYN != 'U' and LorentzianHWHMYN == 'U':
        use_cols = [0,1,2,gamma_hwhm_colid]
        use_names = ['u','l','A','gamma_hwhm']
    else:
        use_cols = [0,1,2]
        use_names = ['u','l','A']
    large_file = is_large_trans_file(trans_filepath)
    trans_reader = read_trans_chunks(trans_filepath, use_cols, use_names)
    desc = 'Processing ' + trans_filename + (' (streaming)' if large_file else '')
    if large_file:
        print('Large transition file detected (>1 GB). Streaming chunks sequentially to reduce memory usage.')
        xsecs = None
        for trans_df_chunk in log_tqdm(trans_reader, desc=desc):
            chunk_xsec = process_exomol_cross_section_chunk(states_part_df,T_list,Tvib_list,Trot_list,P,Q_arr,
                                                     broad,ratio,nbroad,broad_dfs,profile_label,trans_df_chunk,temp_idx)
            xsecs = chunk_xsec if xsecs is None else xsecs + chunk_xsec
        if xsecs is None:
            xsecs = np.zeros_like(wn_grid)
    else:
        trans_chunks = list(trans_reader)
        if len(trans_chunks) == 0:
            xsecs = np.zeros_like(wn_grid)
        else:
            with ProcessPoolExecutor(max_workers=ncputrans) as trans_executor:
                futures = [
                    trans_executor.submit(process_exomol_cross_section_chunk,states_part_df,T_list,Tvib_list,Trot_list,P,Q_arr,
                                          broad,ratio,nbroad,broad_dfs,profile_label,chunk,temp_idx) 
                    for chunk in log_tqdm(trans_chunks, desc=desc)
                ]
                xsecs = np.sum([future.result() for future in log_tqdm(futures, desc='Combining '+trans_filename)], axis=0)
    return xsecs

# Cross sections for ExoMol database
def save_exomol_cross_section(states_part_df, T_list, Tvib_list, Trot_list, P_list, Q_arr):
    """
    Main function to calculate and save cross sections for ExoMol database.

    Processes all transition files, calculates cross sections for each (T, P) combination,
    and saves results in .xsec format files.

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
    """
    from pyexocross.core import (
        read_path,
        profile,
        NLTEMethod,
        abs_emi,
        cutoff,
        UncFilter,
        min_wnl,
        max_wnl,
        data_info,
        wn_grid,
        database,
        ncpufiles,
    )
    print('Calculate cross sections.')
    tot = Timer()
    tot.start()
    broad, ratio, nbroad, broad_dfs = read_broad(read_path)
    # Q = read_exomol_pf(read_path, T)
    profile_label = line_profile(profile)
    
    # Check if profile is pressure-dependent
    pressure_dependent = profile_label not in ['Gaussian', 'Doppler', 'Binned Doppler', 'Binned Gaussion']
    
    # ntemp: L=len(T_list), T/D=len(Trot_list), P=1 (no temperature dimension)
    ntemp = get_ntemp(NLTEMethod, T_list, Trot_list)
    # Determine pressure structure
    # P_list can be: single value, list (one per temp), or list of lists (multiple per temp)
    if isinstance(P_list, (list, np.ndarray)) and len(P_list) > 0:
        first_item = P_list[0]
        if isinstance(first_item, (list, np.ndarray)):
            # Multiple pressures per temperature: P_list is list of lists
            P_per_temp = [list(P_list[i]) if i < len(P_list) else list(P_list[0]) for i in range(ntemp)]
        else:
            # Same pressure list applies to every temperature (cartesian product of T and P)
            P_values = list(P_list)
            P_per_temp = [P_values[:] for _ in range(ntemp)]
    else:
        # Single pressure for all temperatures
        P_single = P_list if isinstance(P_list, (int, float)) else (P_list[0] if len(P_list) > 0 else 1.0)
        P_per_temp = [[P_single] for _ in range(ntemp)]
    
    # Print cross section information once before processing
    print_xsec_info(profile_label, cutoff, UncFilter, min_wnl, max_wnl, 
                    'cm⁻¹', 'cm⁻¹/(molecule cm⁻²)', broad, ratio)
    
    print('Reading transitions and calculating cross sections ...')    
    trans_filepaths = get_part_transfiles(read_path, data_info)
    
    # Process each (T, P) combination separately to save memory
    any_results = False
    xsec_file_count = 0
    for temp_idx in log_tqdm(range(ntemp), desc='\nProcessing cross sections'):
        T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
        
        if pressure_dependent:
            # For pressure-dependent profiles: process each pressure for this temperature
            # This gives num_T * num_P files (for each T, save num_P files)
            for press_idx, P in enumerate(P_per_temp[temp_idx]):
                # Process multiple files in parallel for this (T, P) combination
                with ProcessPoolExecutor(max_workers=ncpufiles) as executor:
                    # Submit reading tasks for each file
                    futures = [executor.submit(process_exomol_cross_section,states_part_df,T_list,Tvib_list,Trot_list,P,Q_arr,broad,ratio,nbroad,broad_dfs,
                                               profile_label,trans_filepath,temp_idx) for trans_filepath in trans_filepaths]
                    xsec = sum([future.result() for future in futures])
                    
                if len(xsec) == 0 or np.all(xsec == 0):
                    print(f'Warning: No cross sections found for T={T} K, P={P} bar. Skipping this combination.')
                    continue
                
                any_results = True
                
                # Save cross sections for this (T, P) combination
                save_xsec_file_plot(wn_grid, xsec, database, profile_label, T, P, temp_idx, Tvib_list, Trot_list)
                xsec_file_count += 1
                print_T_Tvib_Trot_P_path_info(T, Tvib, Trot, P, abs_emi, NLTEMethod, 'Cross sections', None)
                
                # Clear memory
                del xsec
        else:
            # For non-pressure-dependent profiles: process only first pressure (pressure doesn't matter)
            # This gives num_T files (one per temperature)
            P = P_per_temp[temp_idx][0]  # Use first pressure
            
            # Process multiple files in parallel for this temperature
            with ProcessPoolExecutor(max_workers=ncpufiles) as executor:
                # Submit reading tasks for each file
                futures = [executor.submit(process_exomol_cross_section,states_part_df,T_list,Tvib_list,Trot_list,P,Q_arr,broad,ratio,nbroad,broad_dfs,
                                           profile_label,trans_filepath,temp_idx) for trans_filepath in trans_filepaths]
                xsec = sum([future.result() for future in futures])
                
            if len(xsec) == 0 or np.all(xsec == 0):
                print(f'Warning: No cross sections found for T={T} K. Skipping this temperature.')
                continue
            
            any_results = True
            
            # Save cross sections for this temperature
            save_xsec_file_plot(wn_grid, xsec, database, profile_label, T, P, temp_idx, Tvib_list, Trot_list)
            xsec_file_count += 1
            print_T_Tvib_Trot_P_path_info(T, Tvib, Trot, None, abs_emi, NLTEMethod, 'Cross sections', None)
            
            # Clear memory
            del xsec
    
    tot.end()
    print('\nFinished reading transitions and calculating cross sections!\n')
    
    if not any_results:
        raise ValueError("Empty result with the input filter values. Please type new filter values in the input file.")
    
    print(f'All {xsec_file_count} cross sections files have been saved!\n')
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')
