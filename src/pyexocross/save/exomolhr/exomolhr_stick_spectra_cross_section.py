"""
Combined ExoMolHR stick-spectra and cross-section workflow.

Reads and processes the ExoMolHR line list once, then calculates both
stick spectra and cross sections for each temperature/pressure combination.
"""
import numpy as np
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

from pyexocross.base.large_file import save_large_txt
from pyexocross.base.log import (
    log_tqdm,
    print_stick_info,
    print_xsec_info,
    print_file_info,
    print_T_Tvib_Trot_P_path_info,
)
from pyexocross.base.utils import Timer, ensure_dir
from pyexocross.database.load_exomolhr import process_exomolhr_linelist_Q
from pyexocross.plot.plot_stick_spectra import plot_stick_spectra
from pyexocross.plot.plot_cross_section import save_xsec_file_plot
from pyexocross.process.T_n_val import get_ntemp, get_temp_vals
from pyexocross.process.stick_xsec_filepath import stick_spectra_filepath
from pyexocross.save.hitran.hitran_stick_spectra import process_hitran_stick_spectra
from pyexocross.save.exomolhr.exomolhr_cross_section import process_exomolhr_cross_section_chunk
from pyexocross.calculation.calcualte_line_profile import line_profile

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

def process_exomolhr_stick_spectra_cross_section_singleT(A, v, Ep, Epp, gp, T, Q, profile_label, P=None,
                                                         Tvib=None, Trot=None, Evibp=None, Erotp=None, Evibpp=None, Erotpp=None,
                                                         exomolhr_df=None):
    """
    Calculate stick spectra or cross sections for a single temperature/pressure combination.
    """
    if P is None:
        I = process_hitran_stick_spectra(
            A, v, Ep, Epp, gp, T, Q,
            Tvib, Trot, Evibp, Erotp, Evibpp, Erotpp,
        )
        return I, None
    else:
        xsec = process_exomolhr_cross_section_chunk(
            exomolhr_df, T, P, Q, profile_label,
            Tvib, Trot, Evibp, Erotp, Evibpp, Erotpp,
        )
        return None, xsec

def save_exomolhr_stick_spectra_cross_section(exomolhr_df, QNs_col, T_list, Tvib_list, Trot_list, P_list):
    """Run ExoMolHR stick spectra and cross sections using the same line list.

    Reads and processes the data ONCE, then computes both stick spectra and
    cross sections for every temperature (and pressure) combination in parallel.
    """
    from pyexocross.core import (
        NLTEMethod,
        abs_emi,
        QNsformat_list,
        save_path,
        data_info,
        min_wnl,
        max_wnl,
        database,
        PlotStickSpectraYN,
        threshold,
        wn_wl,
        wn_wl_unit,
        UncFilter,
        LTE_NLTE,
        photo,
        wn_grid,
        cutoff,
        profile,
        read_path,
        ncputrans,
    )

    print('Calculate stick spectra and cross sections (read once, calculate both together).\n')
    t_ss = Timer()
    t_xsec = Timer()

    # ---- Prepare data ONCE ----
    print('Preparing ExoMolHR line list data ONCE for all temperatures '
          '(will be used for both stick spectra and cross sections)...')
    A, v, Ep, Epp, gp, Q_arr, Evibp, Erotp, Evibpp, Erotpp = process_exomolhr_linelist_Q(
        exomolhr_df, T_list, Tvib_list, Trot_list,
    )

    profile_label = line_profile(profile)
    pressure_dependent = profile_label not in [
        'Gaussian', 'Doppler', 'Binned Doppler', 'Binned Gaussion',
    ]

    n_temps = get_ntemp(NLTEMethod, T_list, Trot_list)

    # ---- Pressure structure ----
    if isinstance(P_list, (list, np.ndarray)) and len(P_list) > 0:
        first_item = P_list[0]
        if isinstance(first_item, (list, np.ndarray)):
            P_per_temp = [list(P_list[i]) if i < len(P_list) else list(P_list[0])
                          for i in range(n_temps)]
        else:
            P_values = list(P_list)
            P_per_temp = [P_values[:] for _ in range(n_temps)]
    else:
        P_single = (P_list if isinstance(P_list, (int, float))
                     else (P_list[0] if len(P_list) > 0 else 1.0))
        P_per_temp = [[P_single] for _ in range(n_temps)]

    # ---- Absorption / Emission check ----
    if abs_emi == 'Ab':
        print('Absorption stick spectra and cross sections')
    elif abs_emi == 'Em':
        print('Emission stick spectra and cross sections')
    else:
        raise ValueError("Please choose one from: 'Absorption' or 'Emission'.")

    # ---- Stick spectra file format ----
    QNsfmf = (str(QNsformat_list)
              .replace("'", "").replace(",", "").replace("[", "").replace("]", "")
              .replace('d', 's').replace('i', 's').replace('.1f', 's'))
    ss_folder = save_path + 'stick_spectra/files/' + data_info[0] + '/' + database + '/'
    ensure_dir(ss_folder)
    str_min_wnl = str(int(np.floor(min_wnl)))
    str_max_wnl = str(int(np.ceil(max_wnl)))

    if QNsfmf == '':
        ss_fmt = '%15.6f %15.8E %7s %12.6f %7s %12.6f'
    else:
        ss_fmt = '%15.6f %15.8E %7s %12.6f %7s %12.6f ' + QNsfmf + ' ' + QNsfmf

    keep_cols = ['v', 'S', "J'", "E'", 'J"', 'E"'] + QNs_col 
    base_df = exomolhr_df[keep_cols].copy()

    # ---- Read broadening info once (for cross section info line) ----
    from pyexocross.database import read_broad
    broad, ratio, _, _ = read_broad(read_path)

    # ---- Print info ----
    print('Calculate stick spectra and cross sections together.')
    print_stick_info('cm⁻¹', 'cm/molecule')
    print_xsec_info(profile_label, cutoff, UncFilter, min_wnl, max_wnl,
                    'cm⁻¹', 'cm⁻¹/(molecule cm⁻²)', broad, ratio)

    # Prepare tasks for stick spectra and cross sections
    ss_tasks = []
    xsec_tasks = []
    for temp_idx in range(n_temps):
        T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
        Q = Q_arr[temp_idx]
        ss_tasks.append((temp_idx, T, Q, Tvib, Trot))
        
        pressures = P_per_temp[temp_idx] if pressure_dependent else [P_per_temp[temp_idx][0]]
        for P in pressures:
            xsec_tasks.append((temp_idx, T, Q, P, Tvib, Trot))

    # Storage for parallel results
    stick_spectra_results = {}
    
    # ---- Main loop for stick spectra ----
    t_ss.start()
    any_results_ss = False
    ss_file_count = 0
    print('Processing stick spectra in parallel...')
    
    with _executor_context(max_workers=ncputrans) as executor:
        if NLTEMethod in ('L', 'P'):
            futures_ss = {executor.submit(process_exomolhr_stick_spectra_cross_section_singleT,
                                         A, v, Ep, Epp, gp, T, Q, None, None): temp_idx
                         for temp_idx, T, Q, Tvib, Trot in log_tqdm(ss_tasks, desc='\nProcessing stick spectra')}
        else:
            futures_ss = {executor.submit(process_exomolhr_stick_spectra_cross_section_singleT,
                                         A, v, Ep, Epp, gp, T, Q, None, None,
                                         Tvib, Trot, Evibp, Erotp, Evibpp, Erotpp): temp_idx
                         for temp_idx, T, Q, Tvib, Trot in log_tqdm(ss_tasks, desc='\nProcessing stick spectra')}
            
        for future in as_completed(futures_ss):
            temp_idx = futures_ss[future]
            try:
                I, _ = future.result()
                if I is not None and len(I) > 0:
                    stick_spectra_df = base_df.copy()
                    stick_spectra_df['S'] = I
                    stick_spectra_df = stick_spectra_df[['v', 'S', "J'", "E'", 'J"', 'E"'] + QNs_col]
                    
                    if threshold != 'None':
                        stick_spectra_df = stick_spectra_df[stick_spectra_df['S'] >= threshold]
                        
                    if len(stick_spectra_df) > 0:
                        any_results_ss = True
                        stick_spectra_results[temp_idx] = stick_spectra_df
            except Exception as e:
                print(f'Error processing stick spectra for temp_idx={temp_idx}: {e}')

    if any_results_ss:
        print('\nSaving stick spectra results...')
        for temp_idx in list(stick_spectra_results.keys()):
            stick_spectra_df = stick_spectra_results.pop(temp_idx)
            T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
            
            if PlotStickSpectraYN == 'Y':
                plot_stick_spectra(stick_spectra_df, T=T, Tvib=Tvib, Trot=Trot)

            # Unit conversion for saving
            if wn_wl == 'WN':
                unit_fn = 'cm-1__'
            elif wn_wl == 'WL' and wn_wl_unit == 'um':
                stick_spectra_df['v'] = 1e4 / stick_spectra_df['v']
                unit_fn = 'um__'
            elif wn_wl == 'WL' and wn_wl_unit == 'nm':
                stick_spectra_df['v'] = 1e7 / stick_spectra_df['v']
                unit_fn = 'nm__'
            else:
                raise ValueError('Please wirte the unit of wavelength in the input file: um or nm.')
            stick_spectra_df.sort_values(by=['v'], ascending=True, inplace=True)

            ss_path = stick_spectra_filepath(
                ss_folder, T, Tvib, Trot, str_min_wnl, str_max_wnl, unit_fn,
                data_info, wn_wl, UncFilter, threshold, database, abs_emi,
                LTE_NLTE, photo, NLTEMethod,
            )
            save_large_txt(ss_path, stick_spectra_df, fmt=ss_fmt)
            ss_file_count += 1
            print_T_Tvib_Trot_P_path_info(T, Tvib, Trot, None, abs_emi, NLTEMethod,
                                          'Stick spectra', ss_path)
            ss_col = list(stick_spectra_df.columns)
            del stick_spectra_df
            
        print('\nTotal running time for stick spectra:')
        t_ss.end()
        if wn_wl == 'WN':
            ss_col_list = ['Wavenumber'] + ss_col[1:]
            ss_fmt_list = ['%15.6f'] + ss_fmt.split()[1:]
        else:
            ss_col_list = ['Wavelength'] + ss_col[1:]
            ss_fmt_list = ['%15.8E'] + ss_fmt.split()[1:]
        print_file_info('Stick spectra', ss_col_list, ss_fmt_list)
        print(f'\nAll {ss_file_count} stick spectra files have been saved!\n')
        print('* * * * * - - - - - * * * * * - - - * * * * * - - - - - * * * * *\n')
    else:
        print(f'Warning: No transitions found for any temperature. Skipping stick spectra.')

    # ---- Main loop for cross sections ----
    t_xsec.start()
    any_results_xsec = False
    xsec_file_count = 0
    print('Processing cross sections in parallel...')
    
    with _executor_context(max_workers=ncputrans) as executor:
        if NLTEMethod in ('L', 'P'):
            futures_xsec = {executor.submit(process_exomolhr_stick_spectra_cross_section_singleT,
                                           A, v, Ep, Epp, gp, T, Q, profile_label, P,
                                           None, None, None, None, None, None,
                                           exomolhr_df): (temp_idx, P)
                           for temp_idx, T, Q, P, Tvib, Trot in log_tqdm(xsec_tasks, desc='\nProcessing cross sections')}
        else:
            futures_xsec = {executor.submit(process_exomolhr_stick_spectra_cross_section_singleT,
                                           A, v, Ep, Epp, gp, T, Q, profile_label, P,
                                           Tvib, Trot, Evibp, Erotp, Evibpp, Erotpp,
                                           exomolhr_df): (temp_idx, P)
                           for temp_idx, T, Q, P, Tvib, Trot in log_tqdm(xsec_tasks, desc='\nProcessing cross sections')}
            
        for future in as_completed(futures_xsec):
            temp_idx, P = futures_xsec[future]
            T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
            try:
                _, xsec = future.result()
                if xsec is not None and not (len(xsec) == 0 or np.all(xsec == 0)):
                    any_results_xsec = True
                    save_xsec_file_plot(wn_grid, xsec, database, profile_label, T, P,
                                        temp_idx, Tvib_list, Trot_list)
                    xsec_file_count += 1
                    print_T_Tvib_Trot_P_path_info(T, Tvib, Trot,
                                                  P if pressure_dependent else None,
                                                  abs_emi, NLTEMethod, 'Cross sections', None)
                    del xsec
            except Exception as e:
                print(f'Error processing cross sections for T={T} K, P={P} bar: {e}')

    if any_results_xsec:
        print('\nTotal running time for cross sections:')
        t_xsec.end()
        if wn_wl == 'WN':
            print_file_info('Cross sections', ['Wavenumber', 'Cross section'], ['%15.6f', '%15.8E'])
        else:
            print_file_info('Cross sections', ['Wavelength', 'Cross section'], ['%15.8E', '%15.8E'])
        print(f'\nAll {xsec_file_count} cross sections files have been saved!\n')

    if not any_results_ss and not any_results_xsec:
        raise ValueError("Empty result with the input filter values. Please type new filter values in the input file.")

    print('Finished calculating stick spectra and cross sections!\n')
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')
