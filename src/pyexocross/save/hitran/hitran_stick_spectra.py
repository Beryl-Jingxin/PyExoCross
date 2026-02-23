import numpy as np
from concurrent.futures import ThreadPoolExecutor

from pyexocross.base.utils import Timer, ensure_dir
from pyexocross.base.log import log_tqdm, print_stick_info, print_T_Tvib_Trot_P_path_info
from pyexocross.base.large_file import save_large_txt
from pyexocross.calculation.calculate_para import cal_v
from pyexocross.calculation.calculate_intensity import (
    cal_abscoefs,
    cal_abscoefs_nlte_2T
)
from pyexocross.calculation.calculate_emissivity import (
    cal_emicoefs,
    cal_emicoefs_nlte_2T
)
from pyexocross.process.T_n_val import get_ntemp, get_temp_vals
from pyexocross.process.stick_xsec_filepath import stick_spectra_filepath
from pyexocross.database.load_hitran import process_hitran_linelist_Q
from pyexocross.plot.plot_stick_spectra import plot_stick_spectra

def process_hitran_stick_spectra(A, v, Ep, Epp, gp, T, Q, Tvib=None, Trot=None, Evibp=None, Erotp=None, Evibpp=None, Erotpp=None):
    """
    Calculate absorption/emission coefficients for HITRAN stick spectra.

    Unified function that handles both LTE and non-LTE (two-temperature) methods,
    and both absorption and emission calculations.

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
    T : float
        Temperature in Kelvin (for LTE or linewidth calculation)
    Q : float
        Partition function value
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
        Absorption or emission coefficient array, shape (n_lines,)
    """
    from pyexocross.core import abs_emi, NLTEMethod, abundance
    if abs_emi == 'Ab':
        if NLTEMethod == 'L':
            # LTE absorption
            coef = cal_abscoefs([T], [Q], Epp, gp, A, v, abundance)[0, :]
        elif NLTEMethod == 'T':
            # Non-LTE absorption using two temperatures (Tvib, Trot)
            coef = cal_abscoefs_nlte_2T([Tvib], [Trot], [Q], Evibpp, Erotpp, gp, A, v, abundance)[0, :]
        else:
            raise ValueError("Please choose 'LTE' or 'Non-LTE' and 'T' two temperatures (Tvib, Trot) method.")
    elif abs_emi == 'Em':
        if NLTEMethod == 'L':
            # LTE emission
            coef = cal_emicoefs([T], [Q], Ep, gp, A, v, abundance)[0, :]
        elif NLTEMethod == 'T':
            # Non-LTE emission using two temperatures (Tvib, Trot)
            coef = cal_emicoefs_nlte_2T([Tvib], [Trot], [Q], Evibp, Erotp, gp, A, v, abundance)[0, :]
        else:
            raise ValueError("Please choose 'LTE' or 'Non-LTE' and 'T' two temperatures (Tvib, Trot) method.")
    else:
        raise ValueError("Please choose one from: 'Absorption' or 'Emission'.")
    
    return coef

# Stick spectra for HITRAN database
def save_hitran_stick_spectra(hitran_linelist_df, QNs_col, T_list, Tvib_list, Trot_list):
    """
    Main function to calculate and save stick spectra for HITRAN database.

    Processes HITRAN linelist, calculates stick spectra for each temperature,
    and saves results in .stick format files.

    Parameters
    ----------
    hitran_linelist_df : pd.DataFrame
        HITRAN linelist DataFrame with processed quantum numbers
    QNs_col : list of str
        Quantum number column names
    T_list : list of float
        Temperature list for LTE calculations
    Tvib_list : list of float
        Vibrational temperature list for non-LTE two-temperature method (empty [] for LTE)
    Trot_list : list of float
        Rotational temperature list for non-LTE methods (empty [] for LTE)
    """
    from pyexocross.core import (
        abs_emi,
        NLTEMethod,
        QNs_label,
        QNsformat_list,
        QNslabel_list,
        save_path,
        data_info,
        min_wnl,
        max_wnl,
        database,
        PlotStickSpectraYN,
        ncputrans,
        threshold,
        wn_wl,
        wn_wl_unit,
        UncFilter,
        LTE_NLTE,
        photo,
    )

    print('Calculate stick spectra.')  
    print_stick_info('cm⁻¹', 'cm/molecule')
    tot = Timer()
    tot.start()
    
    # Prepare data ONCE using unified function
    print('Preparing HITRAN linelist data ONCE for all temperatures...')
    (A, v, Ep, Epp, gp, n_air, gamma_air, gamma_self, delta_air, 
     Q_arr, Evibp, Erotp, Evibpp, Erotpp) = process_hitran_linelist_Q(hitran_linelist_df, T_list, Tvib_list, Trot_list)
    
    # Print absorption/emission info once
    if abs_emi == 'Ab': 
        print('Absorption stick spectra') 
    elif abs_emi == 'Em': 
        print('Emission stick spectra')
    else:
        raise ValueError("Please choose one from: 'Absorption' or 'Emission'.")
    
    print('Calculating stick spectra ...')    
    
    # Process each temperature separately to save memory
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
    
    # ntemp: L=len(T_list), T/D=len(Trot_list), P=1 (no temperature dimension)
    n_temps = get_ntemp(NLTEMethod, T_list, Trot_list)
    any_results = False
    ss_file_count = 0
    for temp_idx in log_tqdm(range(n_temps), desc='\nProcessing stick spectra'):
        T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
        Q = Q_arr[temp_idx]
        
        # Use threads to keep runtime config available during parallel work.
        with ThreadPoolExecutor(max_workers=ncputrans) as executor:
            future = executor.submit(process_hitran_stick_spectra, 
                                   A, v, Ep, Epp, gp, T, Q, Tvib, Trot, Evibp, Erotp, Evibpp, Erotpp)
            I = future.result()
        
        # Create stick spectra DataFrame with intensity
        stick_spectra_df = base_df.copy()
        stick_spectra_df['S'] = I
        
        # Apply threshold filter if needed
        if threshold != 'None':
            stick_spectra_df = stick_spectra_df[stick_spectra_df['S'] >= threshold]
            
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
        print('\nSaving stick spectra results...')
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
        
    print('\nTotal running time for stick spectra:')
    tot.end()
    print('\nFinished calculating stick spectra!\n')
    
    if not any_results:
        raise ValueError("Empty result with the input filter values. Please type new filter values in the input file.")
    
    print(f'\nAll {ss_file_count} stick spectra files have been saved!\n')
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')
