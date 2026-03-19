"""
Save ExoMolHR stick spectra.
"""
import numpy as np
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

from pyexocross.base.large_file import save_large_txt
from pyexocross.base.log import (
    log_tqdm, 
    print_stick_info, 
    print_file_info, 
    print_T_Tvib_Trot_P_path_info,
)
from pyexocross.base.utils import Timer, ensure_dir
from pyexocross.database.load_exomolhr import process_exomolhr_linelist_Q
from pyexocross.plot.plot_stick_spectra import plot_stick_spectra
from pyexocross.process.T_n_val import get_ntemp, get_temp_vals
from pyexocross.process.stick_xsec_filepath import stick_spectra_filepath
from pyexocross.save.hitran.hitran_stick_spectra import process_hitran_stick_spectra

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

def process_exomolhr_stick_spectra_singleT(A, v, Ep, Epp, gp, T, Q, Tvib=None, Trot=None, Evibp=None, Erotp=None, Evibpp=None, Erotpp=None):
    I = process_hitran_stick_spectra(
        A, v, Ep, Epp, gp, T, Q,
        Tvib, Trot, Evibp, Erotp, Evibpp, Erotpp,
    )
    return I

def save_exomolhr_stick_spectra(exomolhr_df, QNs_col, T_list, Tvib_list, Trot_list):
    """
    Calculate and save stick spectra for ExoMolHR line lists.
    """
    from pyexocross.core import (
        abs_emi,
        NLTEMethod,
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
        ncputrans,
    )

    print('Calculate stick spectra.')
    print_stick_info('cm⁻¹', 'cm/molecule')
    tot = Timer()
    tot.start()

    print('Preparing ExoMolHR line list data ONCE for all temperatures...')
    A, v, Ep, Epp, gp, Q_arr, Evibp, Erotp, Evibpp, Erotpp = process_exomolhr_linelist_Q(
        exomolhr_df,
        T_list,
        Tvib_list,
        Trot_list,
    )

    QNsfmf = (str(QNsformat_list).replace("'","").replace(",","").replace("[","").replace("]","")
              .replace('d','s').replace('i','s').replace('.1f','s'))
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

    n_temps = get_ntemp(NLTEMethod, T_list, Trot_list)
    
    ss_tasks = []
    for temp_idx in range(n_temps):
        T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
        Q = Q_arr[temp_idx]
        ss_tasks.append((temp_idx, T, Q, Tvib, Trot))
        
    stick_spectra_results = {}
    any_results_ss = False
    ss_file_count = 0
    
    print('Processing stick spectra in parallel...')
    with _executor_context(max_workers=ncputrans) as executor:
        if NLTEMethod in ('L', 'P'):
            futures_ss = {executor.submit(process_exomolhr_stick_spectra_singleT,
                                         A, v, Ep, Epp, gp, T, Q, None, None): temp_idx
                         for temp_idx, T, Q, Tvib, Trot in log_tqdm(ss_tasks, desc='\nProcessing stick spectra')}
        else:
            futures_ss = {executor.submit(process_exomolhr_stick_spectra_singleT,
                                         A, v, Ep, Epp, gp, T, Q,
                                         Tvib, Trot, Evibp, Erotp, Evibpp, Erotpp): temp_idx
                         for temp_idx, T, Q, Tvib, Trot in log_tqdm(ss_tasks, desc='\nProcessing stick spectra')}
            
        for future in as_completed(futures_ss):
            temp_idx = futures_ss[future]
            try:
                I = future.result()
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
                ss_folder, T, Tvib, Trot, str_min_wnl, str_max_wnl,
                unit_fn, data_info, wn_wl, UncFilter, threshold,
                database, abs_emi, LTE_NLTE, photo, NLTEMethod,
            )
            save_large_txt(ss_path, stick_spectra_df, fmt=ss_fmt)
            ss_file_count += 1
            print_T_Tvib_Trot_P_path_info(T, Tvib, Trot, None, abs_emi, NLTEMethod, 'Stick spectra', ss_path)
            ss_col = list(stick_spectra_df.columns)
            del stick_spectra_df

    tot.end()
    print('\nFinished calculating stick spectra!\n')

    if not any_results_ss:
        raise ValueError("Empty result with the input filter values. Please type new filter values in the input file.")
        
    if wn_wl == 'WN':
        ss_col_list = ['Wavenumber'] + ss_col[1:]
        ss_fmt_list = ['%15.6f'] + ss_fmt.split()[1:]
    else:
        ss_col_list = ['Wavelength'] + ss_col[1:]
        ss_fmt_list = ['%15.8E'] + ss_fmt.split()[1:]
    print_file_info('Stick spectra', ss_col_list, ss_fmt_list)
    print(f'All {ss_file_count} stick spectra files have been saved!\n')
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')
