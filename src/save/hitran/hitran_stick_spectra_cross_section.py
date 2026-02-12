import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from src.base.utils import Timer, ensure_dir
from src.base.large_file import save_large_txt
from src.base.log import (
    log_tqdm, 
    print_stick_info, 
    print_xsec_info, 
    print_T_Tvib_Trot_P_path_info,
)
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
from src.process.stick_xsec_filepath import stick_spectra_filepath
from src.database.load_hitran import process_hitran_linelist_Q
from src.save.hitran.hitran_stick_spectra import process_hitran_stick_spectra
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
    PseudoLiuLinVoigt,
    PseudoRoccoVoigt
)
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
from src.plot.plot_cross_section import save_xsec_file_plot
from src.plot.plot_stick_spectra import plot_stick_spectra

def process_hitran_stick_spectra_cross_section_singleT(A, v, Ep, Epp, gp, n_air, gamma_air, gamma_self, T, Q, profile_label, P=None,  
                                                       Tvib=None, Trot=None, Evibp=None, Erotp=None, Evibpp=None, Erotpp=None):
    """
    Calculate both stick spectra intensity and cross section for a single temperature.

    Optimized for parallel processing. Calculates absorption/emission coefficients
    and optionally applies line profiles to compute cross sections.

    Parameters
    ----------
    A : np.ndarray
        Einstein A coefficients, shape (n_lines,)
    v : np.ndarray
        Transition wavenumbers, shape (n_lines,)
    Ep : np.ndarray
        Upper state energies, shape (n_lines,)
    Epp : np.ndarray
        Lower state energies, shape (n_lines,)
    gp : np.ndarray
        Upper state degeneracies, shape (n_lines,)
    n_air : np.ndarray
        Temperature exponents, shape (n_lines,)
    gamma_air : np.ndarray
        Air-broadening coefficients, shape (n_lines,)
    gamma_self : np.ndarray
        Self-broadening coefficients, shape (n_lines,)
    T : float
        Temperature in Kelvin
    Q : float
        Partition function value
    profile_label : str or None
        Line profile label (None to skip cross section calculation)
    P : float, optional
        Pressure in bar (required for pressure-dependent profiles)
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
    tuple of (np.ndarray, np.ndarray or None)
        A tuple containing:
        - coef : np.ndarray, absorption/emission coefficient array, shape (n_lines,)
        - xsec : np.ndarray or None, cross-section array if profile_label provided, else None
    """
    from pyexocross.core import threshold, alpha_HWHM, gamma_HWHM, cutoff, wn_grid
    num = len(A)
    if num == 0:
        return None, None
    
    # Calculate absorption/emission coefficients (used for both stick spectra and cross section)
    # Use unified function that handles Ab/Em and LTE/NLTE
    coef = process_hitran_stick_spectra(A, v, Ep, Epp, gp, T, Q, Tvib, Trot, Evibp, Erotp, Evibpp, Erotpp)
    
    # Apply threshold filter using DataFrame (similar to single calculation functions)
    ll_df = pd.DataFrame({'coef': coef, 'v': v, 'n_air': n_air, 'gamma_air': gamma_air, 'gamma_self': gamma_self})
    if threshold != 'None':
        ll_df = ll_df[ll_df['coef'] >= threshold]
        if len(ll_df) == 0:
            return None, None
    coef_filtered = ll_df['coef'].values
    v_filtered = ll_df['v'].values
    n_air_filtered = ll_df['n_air'].values
    gamma_air_filtered = ll_df['gamma_air'].values
    gamma_self_filtered = ll_df['gamma_self'].values
    num_filtered = len(ll_df)
    
    if num_filtered == 0:
        return None, None
    
    # Calculate cross section if profile_label is provided
    xsec = None
    if profile_label is not None:
        if num_filtered > 0:
            # Line profiles
            if profile_label == 'Doppler':
                alpha = DopplerHWHM_alpha(num_filtered, alpha_HWHM, v_filtered, T)
                xsec = cross_section_Doppler(wn_grid, v_filtered, alpha, coef_filtered, cutoff)
            elif profile_label == 'Gaussian':
                alpha = DopplerHWHM_alpha(num_filtered, alpha_HWHM, v_filtered, T)
                xsec = cross_section_Doppler(wn_grid, v_filtered, alpha, coef_filtered, cutoff)
            elif profile_label == 'Lorentzian':
                if P is None:
                    raise ValueError("Pressure is required for Lorentzian profile")
                gamma = LorentzianHWHM_gamma(num_filtered, gamma_HWHM, 0, [], n_air_filtered, gamma_air_filtered, gamma_self_filtered, [], T, P)
                xsec = cross_section_Lorentzian(wn_grid, v_filtered, gamma, coef_filtered, cutoff)
            elif profile_label == 'SciPy Voigt':
                if P is None:
                    raise ValueError("Pressure is required for SciPy Voigt profile")
                alpha = DopplerHWHM_alpha(num_filtered, alpha_HWHM, v_filtered, T)
                sigma = Gaussian_standard_deviation(alpha)    
                gamma = LorentzianHWHM_gamma(num_filtered, gamma_HWHM, 0, [], n_air_filtered, gamma_air_filtered, gamma_self_filtered, [], T, P) 
                xsec = cross_section_SciPyVoigt(wn_grid, v_filtered, sigma, gamma, coef_filtered, cutoff)
            elif profile_label == 'SciPy wofz Voigt':
                if P is None:
                    raise ValueError("Pressure is required for SciPy wofz Voigt profile")
                alpha = DopplerHWHM_alpha(num_filtered, alpha_HWHM, v_filtered, T)
                sigma = Gaussian_standard_deviation(alpha)  
                gamma = LorentzianHWHM_gamma(num_filtered, gamma_HWHM, 0, [], n_air_filtered, gamma_air_filtered, gamma_self_filtered, [], T, P)   
                xsec = cross_section_SciPyWofzVoigt(wn_grid, v_filtered, sigma, gamma, coef_filtered, cutoff)
            elif profile_label == 'Humlicek Voigt':
                if P is None:
                    raise ValueError("Pressure is required for Humlicek Voigt profile")
                alpha = DopplerHWHM_alpha(num_filtered, alpha_HWHM, v_filtered, T)
                gamma = LorentzianHWHM_gamma(num_filtered, gamma_HWHM, 0, [], n_air_filtered, gamma_air_filtered, gamma_self_filtered, [], T, P)
                xsec = cross_section_HumlicekVoigt(wn_grid, v_filtered, alpha, gamma, coef_filtered, cutoff)  
            elif profile_label == 'Thompson pseudo-Voigt':
                if P is None:
                    raise ValueError("Pressure is required for Thompson pseudo-Voigt profile")
                alpha = DopplerHWHM_alpha(num_filtered, alpha_HWHM, v_filtered, T)
                gamma = LorentzianHWHM_gamma(num_filtered, gamma_HWHM, 0, [], n_air_filtered, gamma_air_filtered, gamma_self_filtered, [], T, P)
                eta = PseudoThompsonVoigt(alpha, gamma)
                xsec = cross_section_PseudoVoigt(wn_grid, v_filtered, alpha, gamma, eta, coef_filtered, cutoff)       
            elif profile_label == 'Kielkopf pseudo-Voigt':
                if P is None:
                    raise ValueError("Pressure is required for Kielkopf pseudo-Voigt profile")
                alpha = DopplerHWHM_alpha(num_filtered, alpha_HWHM, v_filtered, T)
                gamma = LorentzianHWHM_gamma(num_filtered, gamma_HWHM, 0, [], n_air_filtered, gamma_air_filtered, gamma_self_filtered, [], T, P)
                eta = PseudoKielkopfVoigt(alpha, gamma)
                xsec = cross_section_PseudoVoigt(wn_grid, v_filtered, alpha, gamma, eta, coef_filtered, cutoff)       
            elif profile_label == 'Olivero pseudo-Voigt':
                if P is None:
                    raise ValueError("Pressure is required for Olivero pseudo-Voigt profile")
                alpha = DopplerHWHM_alpha(num_filtered, alpha_HWHM, v_filtered, T)
                gamma = LorentzianHWHM_gamma(num_filtered, gamma_HWHM, 0, [], n_air_filtered, gamma_air_filtered, gamma_self_filtered, [], T, P)
                eta = PseudoOliveroVoigt(alpha, gamma)
                xsec = cross_section_PseudoVoigt(wn_grid, v_filtered, alpha, gamma, eta, coef_filtered, cutoff)       
            elif profile_label == 'Liu-Lin pseudo-Voigt':
                if P is None:
                    raise ValueError("Pressure is required for Liu-Lin pseudo-Voigt profile")
                alpha = DopplerHWHM_alpha(num_filtered, alpha_HWHM, v_filtered, T)
                gamma = LorentzianHWHM_gamma(num_filtered, gamma_HWHM, 0, [], n_air_filtered, gamma_air_filtered, gamma_self_filtered, [], T, P)
                eta = PseudoLiuLinVoigt(alpha, gamma)
                xsec = cross_section_PseudoVoigt(wn_grid, v_filtered, alpha, gamma, eta, coef_filtered, cutoff)       
            elif profile_label == 'Rocco pseudo-Voigt':
                if P is None:
                    raise ValueError("Pressure is required for Rocco pseudo-Voigt profile")
                alpha = DopplerHWHM_alpha(num_filtered, alpha_HWHM, v_filtered, T)
                gamma = LorentzianHWHM_gamma(num_filtered, gamma_HWHM, 0, [], n_air_filtered, gamma_air_filtered, gamma_self_filtered, [], T, P)
                eta = PseudoRoccoVoigt(alpha, gamma)
                xsec = cross_section_PseudoVoigt(wn_grid, v_filtered, alpha, gamma, eta, coef_filtered, cutoff) 
            elif profile_label == 'Binned Doppler':
                alpha = DopplerHWHM_alpha(num_filtered, alpha_HWHM, v_filtered, T)
                xsec = cross_section_BinnedGaussian(wn_grid, v_filtered, alpha, coef_filtered, cutoff)        
            elif profile_label == 'Binned Gaussion':
                alpha = DopplerHWHM_alpha(num_filtered, alpha_HWHM, v_filtered, T)
                xsec = cross_section_BinnedGaussian(wn_grid, v_filtered, alpha, coef_filtered, cutoff)
            elif profile_label == 'Binned Lorentzian':
                if P is None:
                    raise ValueError("Pressure is required for Binned Lorentzian profile")
                gamma = LorentzianHWHM_gamma(num_filtered, gamma_HWHM, 0, [], n_air_filtered, gamma_air_filtered, gamma_self_filtered, [], T, P)
                xsec = cross_section_BinnedLorentzian(wn_grid, v_filtered, gamma, coef_filtered, cutoff)
            elif profile_label == 'Binned Voigt':
                if P is None:
                    raise ValueError("Pressure is required for Binned Voigt profile")
                alpha = DopplerHWHM_alpha(num_filtered, alpha_HWHM, v_filtered, T)
                sigma = Gaussian_standard_deviation(alpha)     
                gamma = LorentzianHWHM_gamma(num_filtered, gamma_HWHM, 0, [], n_air_filtered, gamma_air_filtered, gamma_self_filtered, [], T, P)
                xsec = cross_section_BinnedVoigt(wn_grid, v_filtered, sigma, gamma, coef_filtered, cutoff)           
            else:
                raise ValueError(f'Unknown profile: {profile_label}')
    
    # Return stick spectra intensity (coef) and cross section
    return coef, xsec

# Stick spectra and cross section for HITRAN database (read data once, calculate both together)
def save_hitran_stick_spectra_cross_section(hitran_linelist_df, QNs_col, T_list, P_list, Tvib_list, Trot_list):
    """
    Combined function to calculate and save both stick spectra and cross sections.

    Reads HITRAN linelist once and processes all temperatures in parallel,
    calculating both stick spectra and cross sections simultaneously.

    Parameters
    ----------
    hitran_linelist_df : pd.DataFrame
        HITRAN linelist DataFrame with processed quantum numbers
    QNs_col : list of str
        Quantum number column names
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
        QNs_label,
        QNsformat_list,
        QNslabel_list,
        save_path,
        data_info,
        min_wnl,
        max_wnl,
        UncFilter,
        threshold,
        cutoff,
        wn_wl,
        wn_wl_unit,
        PlotStickSpectraYN,
        LTE_NLTE,
        photo,
        wn_grid,
        ncputrans,
    )
    print('Calculate stick spectra and cross sections (read once, calculate both together).\n')
    # Separate timers for stick spectra and cross sections
    t_ss = Timer()
    t_xsec = Timer()
    
    # Prepare data ONCE using unified function (extract line data, calculate partition functions, extract Evib/Erot)
    print('Preparing HITRAN linelist data ONCE for all temperatures (will be used for both stick spectra and cross sections)...')
    (A, v, Ep, Epp, gp, n_air, gamma_air, gamma_self, delta_air, 
     Q_arr, Evibp, Erotp, Evibpp, Erotpp) = process_hitran_linelist_Q(hitran_linelist_df, T_list, Tvib_list, Trot_list)
    
    profile_label = line_profile(profile)
    
    # Check if profile is pressure-dependent
    pressure_dependent = profile_label not in ['Gaussian', 'Doppler', 'Binned Doppler', 'Binned Gaussion']
    
    # ntemp: L=len(T_list), T/D=len(Trot_list), P=1 (no temperature dimension)
    ntemp = get_ntemp(NLTEMethod, T_list, Trot_list)
    # Determine pressure structure
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
    
    # Print absorption/emission info once
    if abs_emi == 'Ab': 
        print('Absorption stick spectra and cross sections') 
    elif abs_emi == 'Em': 
        print('Emission stick spectra and cross sections')
    else:
        raise ValueError("Please choose one from: 'Absorption' or 'Emission'.")
    
    # Prepare file format strings for stick spectra
    QN_format_noJ = [QNsformat_list[i] for i in [QNslabel_list.index(j) for j in QNs_label]]
    QNsfmf = (str(QN_format_noJ).replace("'","").replace(",","").replace("[","").replace("]","")
              .replace('d','s').replace('i','s').replace('.1f','s'))
    ss_folder = save_path + 'stick_spectra/files/'+data_info[0]+'/'+database+'/'
    ensure_dir(ss_folder)
    str_min_wnl = str(int(np.floor(min_wnl)))
    str_max_wnl = str(int(np.ceil(max_wnl)))
    
    if QNsfmf == '':
        ss_fmt = '%12.8E %12.8E %7s %12.4f %7s %12.4f'
    else:
        ss_fmt = '%12.8E %12.8E %7s %12.4f %7s %12.4f ' + QNsfmf + ' ' + QNsfmf
    
    # Prepare base DataFrame columns (extract once)
    ss_colname = ['v','S',"J'","E'",'J"','E"'] + QNs_col
    base_df = hitran_linelist_df[ss_colname].copy()
    
    # Initialize results storage
    # Stick spectra: once per temperature (pressure doesn't matter)
    stick_spectra_results = {}
    # Cross sections: for pressure-dependent profiles, store per (T,P) combination
    # For non-pressure-dependent, store once per temperature
    if pressure_dependent:
        xsec_results = {}  # Store as dict: (temp_idx, P) -> xsec array
    else:
        xsec_results = {}  # Store as dict: temp_idx -> xsec array
    
    # Print info
    print('Calculate stick spectra and cross sections together.')
    print_stick_info('cm⁻¹', 'cm/molecule')
    print_xsec_info(profile_label, cutoff, UncFilter, min_wnl, max_wnl, 
                    'cm⁻¹', 'cm⁻¹/(molecule cm⁻²)', [], [])
    
    # Process all temperatures in a single loop, calculating both stick spectra and cross sections
    # Use parallel processing to speed up calculations
    t_ss.start()
    t_xsec.start()
    any_results_ss = False
    any_results_xsec = False
    ss_file_count = 0
    xsec_file_count = 0
    
    # Prepare tasks for parallel processing
    # For stick spectra: one task per temperature (ntemp: L/T/D/P)
    ss_tasks = []
    for temp_idx in range(ntemp):
        T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
        Q = Q_arr[temp_idx]
        ss_tasks.append((temp_idx, T, Q, Tvib, Trot))
    
    # For cross sections: prepare tasks based on pressure dependency
    if pressure_dependent:
        # For pressure-dependent profiles: one task per (T, P) combination
        xsec_tasks = []
        for temp_idx in range(ntemp):
            T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
            Q = Q_arr[temp_idx]
            P_list_for_temp = P_per_temp[temp_idx] if temp_idx < len(P_per_temp) else P_per_temp[0] if len(P_per_temp) > 0 else [1.0]
            for P in P_list_for_temp:
                xsec_tasks.append((temp_idx, T, Q, P, Tvib, Trot))
    else:
        # For non-pressure-dependent profiles: one task per temperature
        xsec_tasks = []
        for temp_idx in range(ntemp):
            T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
            Q = Q_arr[temp_idx]
            P = P_per_temp[temp_idx][0] if temp_idx < len(P_per_temp) and len(P_per_temp[temp_idx]) > 0 else (P_per_temp[0][0] if len(P_per_temp) > 0 and len(P_per_temp[0]) > 0 else 1.0)
            xsec_tasks.append((temp_idx, T, Q, P, Tvib, Trot))

    # Process stick spectra in parallel
    print('Processing stick spectra in parallel...')
    with ProcessPoolExecutor(max_workers=ncputrans) as executor:
        if NLTEMethod == 'L' or NLTEMethod == 'P':
            futures_ss = {executor.submit(process_hitran_stick_spectra_cross_section_singleT,
                                          A, v, Ep, Epp, gp, n_air, gamma_air, gamma_self, T, Q, None, None): temp_idx
                         for temp_idx, T, Q, Tvib, Trot in log_tqdm(ss_tasks, desc='\nProcessing stick spectra')}
        else:
            futures_ss = {executor.submit(process_hitran_stick_spectra_cross_section_singleT,
                                          A, v, Ep, Epp, gp, n_air, gamma_air, gamma_self, T, Q, None, None,
                                          Tvib, Trot, Evibp, Erotp, Evibpp, Erotpp): temp_idx
                         for temp_idx, T, Q, Tvib, Trot in log_tqdm(ss_tasks, desc='\nProcessing stick spectra')}
        
        for future in as_completed(futures_ss):
            temp_idx = futures_ss[future]
            try:
                I, _ = future.result()
                
                if I is not None and len(I) > 0:
                    # Create stick spectra DataFrame with intensity
                    stick_spectra_df = base_df.copy()
                    stick_spectra_df['S'] = I
                    
                    # Apply threshold filter if needed
                    if threshold != 'None':
                        stick_spectra_df = stick_spectra_df[stick_spectra_df['S'] >= threshold]
                        
                    if len(stick_spectra_df) > 0:
                        any_results_ss = True
                        stick_spectra_results[temp_idx] = stick_spectra_df
            except Exception as e:
                print(f'Error processing stick spectra for temp_idx={temp_idx}: {e}')
    
    # Process cross sections in parallel
    print('Processing cross sections in parallel...')
    with ProcessPoolExecutor(max_workers=ncputrans) as executor:
        if NLTEMethod in ('L', 'P'):
            if pressure_dependent:
                futures_xsec = {executor.submit(process_hitran_stick_spectra_cross_section_singleT,
                                               A, v, Ep, Epp, gp, n_air, gamma_air, gamma_self, T, Q, profile_label, P): (temp_idx, P)
                               for temp_idx, T, Q, P, Tvib, Trot in log_tqdm(xsec_tasks, desc='\nProcessing cross sections')}
            else:
                futures_xsec = {executor.submit(process_hitran_stick_spectra_cross_section_singleT,
                                               A, v, Ep, Epp, gp, n_air, gamma_air, gamma_self, T, Q, profile_label, P): temp_idx
                               for temp_idx, T, Q, P, Tvib, Trot in log_tqdm(xsec_tasks, desc='\nProcessing cross sections')}
        elif NLTEMethod in ('T', 'D'):
            if pressure_dependent:
                futures_xsec = {executor.submit(process_hitran_stick_spectra_cross_section_singleT,
                                               A, v, Ep, Epp, gp, n_air, gamma_air, gamma_self, T, Q, profile_label, P,
                                               Tvib, Trot, Evibp, Erotp, Evibpp, Erotpp): (temp_idx, P)
                              for temp_idx, T, Q, P, Tvib, Trot in log_tqdm(xsec_tasks, desc='\nProcessing cross sections')}
            else:
                futures_xsec = {executor.submit(process_hitran_stick_spectra_cross_section_singleT,
                                               A, v, Ep, Epp, gp, n_air, gamma_air, gamma_self, T, Q, profile_label, P,
                                               Tvib, Trot, Evibp, Erotp, Evibpp, Erotpp): temp_idx
                              for temp_idx, T, Q, P, Tvib, Trot in log_tqdm(xsec_tasks, desc='\nProcessing cross sections')}
        else:
            raise ValueError("Please choose: 'LTE' or 'Non-LTE' (T, D or P) method.")
        
        for future in as_completed(futures_xsec):
            try:
                _, xsec = future.result()
                
                if pressure_dependent:
                    temp_idx, P = futures_xsec[future]
                    if xsec is not None and len(xsec) > 0 and not np.all(xsec == 0):
                        any_results_xsec = True
                        xsec_results[(temp_idx, P)] = xsec
                else:
                    temp_idx = futures_xsec[future]
                    if xsec is not None and len(xsec) > 0 and not np.all(xsec == 0):
                        any_results_xsec = True
                        xsec_results[temp_idx] = xsec
            except Exception as e:
                if pressure_dependent:
                    temp_idx, P = futures_xsec[future]
                    print(f'Error processing cross sections for temp_idx={temp_idx}, P={P}: {e}')
                else:
                    temp_idx = futures_xsec[future]
                    print(f'Error processing cross sections for temp_idx={temp_idx}: {e}')
    
    # Save stick spectra results
    if any_results_ss:
        print('\nSaving stick spectra results...')
        for temp_idx, stick_spectra_df in stick_spectra_results.items():
            T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
            
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
                plot_stick_spectra(stick_spectra_df_plot, T=T, Tvib=Tvib, Trot=Trot)
                del stick_spectra_df_plot
        
        print('\nTotal running time for stick spectra:')
        t_ss.end()
        print(f'\nAll {ss_file_count} stick spectra files have been saved!\n')
        print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')
    
    # Save cross section results
    if any_results_xsec:
        print('\nSaving cross sections results...')
        if pressure_dependent:
            # Save for each (T, P) combination
            for (temp_idx, P), xsec in xsec_results.items():
                T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
                save_xsec_file_plot(wn_grid, xsec, database, profile_label, T, P, temp_idx, Tvib_list, Trot_list)
                xsec_file_count += 1
                print_T_Tvib_Trot_P_path_info(T, Tvib, Trot, P, abs_emi, NLTEMethod, 'Cross sections', None)
                del xsec
        else:
            # Save for each temperature
            for temp_idx, xsec in xsec_results.items():
                T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
                P = P_per_temp[temp_idx][0]
                save_xsec_file_plot(wn_grid, xsec, database, profile_label, T, P, temp_idx, Tvib_list, Trot_list)
                xsec_file_count += 1
                print_T_Tvib_Trot_P_path_info(T, Tvib, Trot, None, abs_emi, NLTEMethod, 'Cross sections', None)
                del xsec
        
        print('\nTotal running time for cross sections:')
        t_xsec.end()
    
    # Check if we have any results
    if not any_results_ss and not any_results_xsec:
        raise ValueError("Empty result with the input filter values. Please type new filter values in the input file.")
    
    if any_results_ss:
        print(f'\nAll {ss_file_count} stick spectra files have been saved!\n')
    if any_results_xsec:
        print(f'All {xsec_file_count} cross sections files have been saved!\n')
    print('Finished calculating stick spectra and cross sections!\n')
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')
