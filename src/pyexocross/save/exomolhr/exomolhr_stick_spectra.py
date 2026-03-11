"""
Save ExoMolHR stick spectra.
"""
import numpy as np

from pyexocross.base.large_file import save_large_txt
from pyexocross.base.log import log_tqdm, print_stick_info, print_T_Tvib_Trot_P_path_info
from pyexocross.base.utils import Timer, ensure_dir
from pyexocross.database.load_exomolhr import process_exomolhr_linelist_Q
from pyexocross.plot.plot_stick_spectra import plot_stick_spectra
from pyexocross.process.T_n_val import get_ntemp, get_temp_vals
from pyexocross.process.stick_xsec_filepath import stick_spectra_filepath
from pyexocross.save.hitran.hitran_stick_spectra import process_hitran_stick_spectra


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

    QNsfmf = (str(QNsformat_list + QNsformat_list).replace("'","").replace(",","").replace("[","").replace("]","")
              .replace('d','s').replace('i','s').replace('.1f','s'))
    ss_folder = save_path + 'stick_spectra/files/' + data_info[0] + '/' + database + '/'
    ensure_dir(ss_folder)
    str_min_wnl = str(int(np.floor(min_wnl)))
    str_max_wnl = str(int(np.ceil(max_wnl)))

    if QNsfmf == '':
        ss_fmt = '%12.8E %12.8E %7s %12.4f %7s %12.4f'
    else:
        ss_fmt = '%12.8E %12.8E %7s %12.4f %7s %12.4f ' + QNsfmf

    keep_cols = ['v', "J'", "E'", 'J"', 'E"'] + QNs_col
    base_df = exomolhr_df[keep_cols].copy()

    n_temps = get_ntemp(NLTEMethod, T_list, Trot_list)
    any_results = False
    ss_file_count = 0
    for temp_idx in log_tqdm(range(n_temps), desc='\nProcessing stick spectra'):
        T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
        Q = Q_arr[temp_idx]
        I = process_hitran_stick_spectra(
            A,
            v,
            Ep,
            Epp,
            gp,
            T,
            Q,
            Tvib,
            Trot,
            Evibp,
            Erotp,
            Evibpp,
            Erotpp,
        )

        stick_spectra_df = base_df.copy()
        stick_spectra_df['S'] = I
        stick_spectra_df = stick_spectra_df[['v', 'S', "J'", "E'", 'J"', 'E"'] + QNs_col]
        if threshold != 'None':
            stick_spectra_df = stick_spectra_df[stick_spectra_df['S'] >= threshold]

        if len(stick_spectra_df) == 0:
            print(f'Warning: No transitions found for T={T} K. Skipping this temperature.')
            continue

        any_results = True

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
            ss_folder,
            T,
            Tvib,
            Trot,
            str_min_wnl,
            str_max_wnl,
            unit_fn,
            data_info,
            wn_wl,
            UncFilter,
            threshold,
            database,
            abs_emi,
            LTE_NLTE,
            photo,
            NLTEMethod,
        )
        save_large_txt(ss_path, stick_spectra_df, fmt=ss_fmt)
        ss_file_count += 1
        print_T_Tvib_Trot_P_path_info(T, Tvib, Trot, None, abs_emi, NLTEMethod, 'Stick spectra', ss_path)

        if PlotStickSpectraYN == 'Y':
            stick_spectra_df_plot = stick_spectra_df.copy()
            if wn_wl == 'WL':
                if wn_wl_unit == 'um':
                    stick_spectra_df_plot['v'] = 1e4 / stick_spectra_df_plot['v']
                elif wn_wl_unit == 'nm':
                    stick_spectra_df_plot['v'] = 1e7 / stick_spectra_df_plot['v']
            plot_stick_spectra(stick_spectra_df_plot, T=T, Tvib=Tvib, Trot=Trot)

    tot.end()
    print('\nFinished calculating stick spectra!\n')

    if not any_results:
        raise ValueError("Empty result with the input filter values. Please type new filter values in the input file.")

    print(f'All {ss_file_count} stick spectra files have been saved!\n')
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')
