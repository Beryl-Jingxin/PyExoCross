import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

from src.base.utils import Timer
from src.base.log import log_tqdm, print_xsec_info, print_T_Tvib_Trot_P_path_info
from base.config_manager import get_config
(
    DopplerHWHMYN,
    LorentzianHWHMYN,
    database,
    predissocYN,
    check_predissoc,
    profile,
) = get_config()
from src.process.T_n_val import get_ntemp, get_temp_vals
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
from src.calculation.calcualte_line_profile import (
    DopplerHWHM_alpha,
    Gaussian_standard_deviation,
    LorentzianHWHM_gamma,
    line_profile,
    PseudoThompsonVoigt,
    PseudoKielkopfVoigt,
    PseudoOliveroVoigt,
    PseudoLiuLinVoigt,
    PseudoRoccoVoigt,
)
from src.database.load_hitran import process_hitran_linelist, process_hitran_linelist_Q
from src.save.hitran.hitran_stick_spectra import process_hitran_stick_spectra
from src.plot.plot_cross_section import save_xsec_file_plot




def process_hitran_cross_section_chunk(hitran_linelist_df, T, P, Q, profile_label, Tvib=None, Trot=None, Evibp=None, Erotp=None, Evibpp=None, Erotpp=None):
    """
    Calculate cross section for a single temperature (optimized for parallel processing).

    Processes HITRAN linelist, calculates line strengths, applies line profiles,
    and computes cross sections on the wavenumber grid.

    Parameters
    ----------
    hitran_linelist_df : pd.DataFrame
        HITRAN linelist DataFrame with processed quantum numbers
    T : float
        Temperature in Kelvin
    P : float
        Pressure in bar
    Q : float
        Partition function value
    profile_label : str
        Line profile label (e.g., 'Voigt', 'Doppler', 'Lorentzian')
    Tvib : float, optional
        Vibrational temperature in Kelvin (for non-LTE T mode only)
    Trot : float, optional
        Rotational temperature in Kelvin (for non-LTE T mode only)
    Evibp : np.ndarray, optional
        Upper state vibrational energies (for non-LTE T mode only)
    Erotp : np.ndarray, optional
        Upper state rotational energies (for non-LTE T mode only)
    Evibpp : np.ndarray, optional
        Lower state vibrational energies (for non-LTE T mode only)
    Erotpp : np.ndarray, optional
        Lower state rotational energies (for non-LTE T mode only)

    Returns
    -------
    np.ndarray
        Cross-section array, shape (n_points,)
    """
    from pyexocross.core import threshold, alpha_HWHM, gamma_HWHM, wn_grid, cutoff

    A, v, Ep, Epp, gp, n_air, gamma_air, gamma_self, delta_air = process_hitran_linelist(hitran_linelist_df)
    num = len(A)
    if num > 0:
        # Use unified function that handles Ab/Em and LTE/NLTE
        coef = process_hitran_stick_spectra(A, v, Ep, Epp, gp, T, Q, Tvib, Trot, Evibp, Erotp, Evibpp, Erotpp)  
        ll_df = pd.DataFrame({'coef':coef, 'v':v, 'n_air':n_air, 'gamma_air':gamma_air, 'gamma_self':gamma_self})
        if threshold != 'None':
            ll_df = ll_df[ll_df['coef'] >= threshold]  
            coef = ll_df['coef'].values
            v = ll_df['v'].values
            n_air = ll_df['n_air'].values
            gamma_air = ll_df['gamma_air'].values
            gamma_self = ll_df['gamma_self'].values
            num = len(ll_df)         
        else:
            pass   
        if num > 0:
            # Line profiles
            if profile_label == 'Doppler':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                xsec = cross_section_Doppler(wn_grid, v, alpha, coef, cutoff)
            elif profile_label == 'Gaussian':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                xsec = cross_section_Doppler(wn_grid, v, alpha, coef, cutoff)
            elif profile_label == 'Lorentzian':
                gamma = LorentzianHWHM_gamma(num, gamma_HWHM, 0, [], n_air, gamma_air, gamma_self, [], T, P)
                xsec = cross_section_Lorentzian(wn_grid, v, gamma, coef, cutoff)
            elif profile_label == 'SciPy Voigt':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                sigma = Gaussian_standard_deviation(alpha)    
                gamma = LorentzianHWHM_gamma(num, gamma_HWHM, 0, [], n_air, gamma_air, gamma_self, [], T, P) 
                xsec = cross_section_SciPyVoigt(wn_grid, v, sigma, gamma, coef, cutoff)
            elif profile_label == 'SciPy wofz Voigt':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                sigma = Gaussian_standard_deviation(alpha)  
                gamma = LorentzianHWHM_gamma(num, gamma_HWHM, 0, [], n_air, gamma_air, gamma_self, [], T, P)   
                xsec = cross_section_SciPyWofzVoigt(wn_grid, v, sigma, gamma, coef, cutoff)
            elif profile_label == 'Humlicek Voigt':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                gamma = LorentzianHWHM_gamma(num, gamma_HWHM, 0, [], n_air, gamma_air, gamma_self, [], T, P)
                xsec = cross_section_HumlicekVoigt(wn_grid, v, alpha, gamma, coef, cutoff)  
            elif profile_label == 'Thompson pseudo-Voigt':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                gamma = LorentzianHWHM_gamma(num, gamma_HWHM, 0, [], n_air, gamma_air, gamma_self, [], T, P)
                eta = PseudoThompsonVoigt(alpha, gamma)
                xsec = cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, coef, cutoff)       
            elif profile_label == 'Kielkopf pseudo-Voigt':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                gamma = LorentzianHWHM_gamma(num, gamma_HWHM, 0, [], n_air, gamma_air, gamma_self, [], T, P)
                eta = PseudoKielkopfVoigt(alpha, gamma)
                xsec = cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, coef, cutoff)       
            elif profile_label == 'Olivero pseudo-Voigt':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                gamma = LorentzianHWHM_gamma(num, gamma_HWHM, 0, [], n_air, gamma_air, gamma_self, [], T, P)
                eta = PseudoOliveroVoigt(alpha, gamma)
                xsec = cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, coef, cutoff)       
            elif profile_label == 'Liu-Lin pseudo-Voigt':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                gamma = LorentzianHWHM_gamma(num, gamma_HWHM, 0, [], n_air, gamma_air, gamma_self, [], T, P)
                eta = PseudoLiuLinVoigt(alpha, gamma)
                xsec = cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, coef, cutoff)       
            elif profile_label == 'Rocco pseudo-Voigt':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                gamma = LorentzianHWHM_gamma(num, gamma_HWHM, 0, [], n_air, gamma_air, gamma_self, [], T, P)
                eta = PseudoRoccoVoigt(alpha, gamma)
                xsec = cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, coef, cutoff) 
            elif profile_label == 'Binned Doppler':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                xsec = cross_section_BinnedGaussian(wn_grid, v, alpha, coef, cutoff)        
            elif profile_label == 'Binned Gaussion':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                xsec = cross_section_BinnedGaussian(wn_grid, v, alpha, coef, cutoff)
            elif profile_label == 'Binned Lorentzian':
                gamma = LorentzianHWHM_gamma(num, gamma_HWHM, 0, [], n_air, gamma_air, gamma_self, [], T, P)
                xsec = cross_section_BinnedLorentzian(wn_grid, v, gamma, coef, cutoff)
            elif profile_label == 'Binned Voigt':
                alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
                sigma = Gaussian_standard_deviation(alpha)     
                gamma = LorentzianHWHM_gamma(num, gamma_HWHM, 0, [], n_air, gamma_air, gamma_self, [], T, P)
                xsec = cross_section_BinnedVoigt(wn_grid, v, sigma, gamma, coef, cutoff)           
            else:
                raise ValueError('Please choose line profile from the list.')  
        else:      
                xsec = np.zeros_like(wn_grid)
    else:
        xsec = np.zeros_like(wn_grid)
    return xsec

# Cross sections for HITRAN database
def save_hitran_cross_section(hitran_linelist_df, T_list, P_list, Tvib_list, Trot_list):
    """
    Main function to calculate and save cross sections for HITRAN database.

    Processes HITRAN linelist, calculates cross sections for each (T, P) combination,
    and saves results in .xsec format files.

    Parameters
    ----------
    hitran_linelist_df : pd.DataFrame
        HITRAN linelist DataFrame with processed quantum numbers
    T_list : list of float
        Temperature list for LTE calculations
    P_list : float, list, or list of lists
        Pressure(s) - single value, flat list, or list of lists per temperature
    Tvib_list : list of float
        Vibrational temperature list for non-LTE two-temperature method (empty [] for LTE)
    Trot_list : list of float
        Rotational temperature list for non-LTE methods (empty [] for LTE)
    """
    from pyexocross.core import (
        NLTEMethod,
        abs_emi,
        cutoff,
        UncFilter,
        min_wnl,
        max_wnl,
        wn_grid,
        ncputrans,
    )

    print('Calculating cross sections ...')  
    t = Timer()
    t.start()

    # Prepare data ONCE using unified function
    print('Preparing HITRAN linelist data ONCE for all temperatures...')
    (A, v, Ep, Epp, gp, n_air, gamma_air, gamma_self, delta_air, 
     Q_arr, Evibp, Erotp, Evibpp, Erotpp) = process_hitran_linelist_Q(hitran_linelist_df, T_list, Tvib_list, Trot_list)
    profile_label = line_profile(profile)
    print(profile_label, 'profile')
    
    # Check if profile is pressure-dependent
    pressure_dependent = profile_label not in ['Gaussian', 'Doppler', 'Binned Doppler', 'Binned Gaussion']
    
    # ntemp: L=len(T_list), T/D=len(Trot_list), P=1 (no temperature dimension)
    n_temps = get_ntemp(NLTEMethod, T_list, Trot_list)
    # Determine pressure structure
    # P_list can be: single value, list (one per temp), or list of lists (multiple per temp)
    if isinstance(P_list, (list, np.ndarray)) and len(P_list) > 0:
        first_item = P_list[0]
        if isinstance(first_item, (list, np.ndarray)):
            # Multiple pressures per temperature: P_list is list of lists
            P_per_temp = [list(P_list[i]) if i < len(P_list) else list(P_list[0]) for i in range(n_temps)]
        else:
            # Same pressure list applies to every temperature (cartesian product of T and P)
            P_values = list(P_list)
            P_per_temp = [P_values[:] for _ in range(n_temps)]
    else:
        # Single pressure for all temperatures
        P_single = P_list if isinstance(P_list, (int, float)) else (P_list[0] if len(P_list) > 0 else 1.0)
        P_per_temp = [[P_single] for _ in range(n_temps)]
    
    # Print cross section information once before processing
    print_xsec_info(profile_label, cutoff, UncFilter, min_wnl, max_wnl, 
                    'cm⁻¹', 'cm⁻¹/(molecule cm⁻²)', [], [])
    
    print('Calculating cross sections ...')    
    
    # Process each (T, P) combination separately to save memory
    any_results = False
    xsec_file_count = 0
    temp_info = [
        (temp_idx, *get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list), Q_arr[temp_idx])
        for temp_idx in range(n_temps)
    ]
    temp_info_map = {temp_idx: (T, Tvib, Trot, Q) for temp_idx, T, Tvib, Trot, Q in temp_info}
    
    if pressure_dependent:
        # For pressure-dependent profiles: process each (T, P) combination
        # This gives num_T * num_P files (for each T, save num_P files)
        tp_combinations = [
            (temp_idx, T, P, Q, Tvib, Trot)
            for temp_idx, T, Tvib, Trot, Q in temp_info
            for P in P_per_temp[temp_idx]
        ]
        
        # Process all (T, P) combinations in parallel
        with ProcessPoolExecutor(max_workers=ncputrans) as executor:
            if NLTEMethod in ('L', 'P'):
                futures = {
                    executor.submit(
                        process_hitran_cross_section_chunk,
                        hitran_linelist_df, T, P, Q, profile_label
                    ): (temp_idx, T, P)
                    for temp_idx, T, P, Q, Tvib, Trot in log_tqdm(tp_combinations, desc='\nProcessing cross sections')
                }
            elif NLTEMethod in ('T', 'D'):
                futures = {
                    executor.submit(
                        process_hitran_cross_section_chunk,
                        hitran_linelist_df, T, P, Q, profile_label,
                        Tvib, Trot, Evibp, Erotp, Evibpp, Erotpp
                    ): (temp_idx, T, P)
                    for temp_idx, T, P, Q, Tvib, Trot in log_tqdm(tp_combinations, desc='\nProcessing cross sections')
                }
            else:
                raise ValueError("Please choose: 'LTE' or 'Non-LTE' (T, D or P) method.")
            
            for future in as_completed(futures):
                temp_idx, T, P = futures[future]
                T, Tvib, Trot, _Q = temp_info_map[temp_idx]
                try:
                    xsec = future.result()
                    
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
                except Exception as e:
                    print(f'Error processing T={T} K, P={P} bar: {e}')
    else:
        # For non-pressure-dependent profiles: process only first pressure (pressure doesn't matter)
        # This gives num_T files (one per temperature)
        # Process all temperatures in parallel
        with ProcessPoolExecutor(max_workers=ncputrans) as executor:
            xsec_tasks = [
                (temp_idx, T, Q, P_per_temp[temp_idx][0], Tvib, Trot)
                for temp_idx, T, Tvib, Trot, Q in temp_info
            ]
            if NLTEMethod in ('L', 'P'):
                futures = {
                    executor.submit(
                        process_hitran_cross_section_chunk,
                        hitran_linelist_df, T, P, Q, profile_label
                    ): temp_idx
                    for temp_idx, T, Q, P, Tvib, Trot in log_tqdm(xsec_tasks, desc='\nProcessing cross sections')
                }
            else:
                futures = {
                    executor.submit(
                        process_hitran_cross_section_chunk,
                        hitran_linelist_df, T, P, Q, profile_label,
                        Tvib, Trot, Evibp, Erotp, Evibpp, Erotpp
                    ): temp_idx
                    for temp_idx, T, Q, P, Tvib, Trot in log_tqdm(xsec_tasks, desc='\nProcessing cross sections')
                }

            for future in as_completed(futures):
                temp_idx = futures[future]
                T, Tvib, Trot, _Q = temp_info_map[temp_idx]
                P = P_per_temp[temp_idx][0]
                try:
                    xsec = future.result()
                    
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
                except Exception as e:
                    print(f'Error processing T={T} K: {e}')
    
    t.end()  
    print('\nFinished calculating cross sections!\n')
    
    if not any_results:
        raise ValueError("Empty result with the input filter values. Please type new filter values in the input file.")
    
    print(f'All {xsec_file_count} cross sections files have been saved!\n')
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')
