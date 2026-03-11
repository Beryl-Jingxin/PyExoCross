"""
Save ExoMolHR cross sections.
"""
import numpy as np

from pyexocross.base.log import log_tqdm, print_T_Tvib_Trot_P_path_info, print_xsec_info
from pyexocross.base.utils import Timer
from pyexocross.calculation.calculate_cross_section import (
    cross_section_BinnedGaussian,
    cross_section_BinnedLorentzian,
    cross_section_BinnedVoigt,
    cross_section_Doppler,
    cross_section_HumlicekVoigt,
    cross_section_Lorentzian,
    cross_section_PseudoVoigt,
    cross_section_SciPyVoigt,
    cross_section_SciPyWofzVoigt,
)
from pyexocross.calculation.calcualte_line_profile import (
    DopplerHWHM_alpha,
    Gaussian_standard_deviation,
    LorentzianHWHM_gamma,
    PseudoKielkopfVoigt,
    PseudoLiuLinVoigt,
    PseudoOliveroVoigt,
    PseudoRoccoVoigt,
    PseudoThompsonVoigt,
    line_profile,
)
from pyexocross.database import extract_broad, read_broad
from pyexocross.database.load_exomolhr import process_exomolhr_linelist_Q
from pyexocross.plot.plot_cross_section import save_xsec_file_plot
from pyexocross.process.T_n_val import get_ntemp, get_temp_vals
from pyexocross.save.hitran.hitran_stick_spectra import process_hitran_stick_spectra


def process_exomolhr_cross_section_chunk(exomolhr_df, T, P, Q, profile_label,
                                         Tvib=None, Trot=None,
                                         Evibp=None, Erotp=None, Evibpp=None, Erotpp=None):
    """Calculate one ExoMolHR cross section."""
    from pyexocross.core import threshold, alpha_HWHM, gamma_HWHM, wn_grid, cutoff, read_path

    A = exomolhr_df['A'].values
    v = exomolhr_df['v'].values
    Ep = exomolhr_df["E'"].values
    Epp = exomolhr_df['E"'].values
    gp = exomolhr_df["g'"].values
    num = len(A)
    if num == 0:
        return np.zeros_like(wn_grid)

    coef = process_hitran_stick_spectra(
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
    ll_df = exomolhr_df[['v', 'J"']].copy()
    ll_df['coef'] = coef
    if threshold != 'None':
        ll_df = ll_df[ll_df['coef'] >= threshold]
    if len(ll_df) == 0:
        return np.zeros_like(wn_grid)

    v = ll_df['v'].values
    coef = ll_df['coef'].values
    num = len(ll_df)

    if 'DOP' not in profile_label.upper() and 'GAU' not in profile_label.upper():
        broad, ratio, nbroad, broad_dfs = read_broad(read_path)
        gamma_L = {}
        n_air = {}
        for i in range(nbroad):
            if broad[i] == 'Default':
                gamma_L[i] = np.full(num, float(broad_dfs[i]['gamma_L'][0]), dtype=np.float64) * ratio[i]
                n_air[i] = np.full(num, float(broad_dfs[i]['n_air'][0]), dtype=np.float64) * ratio[i]
            else:
                broad_gamma, broad_nair = extract_broad(broad_dfs[i], ll_df)
                gamma_L[i] = np.asarray(broad_gamma, dtype=np.float64) * ratio[i]
                n_air[i] = np.asarray(broad_nair, dtype=np.float64) * ratio[i]
        gamma = LorentzianHWHM_gamma(num, gamma_HWHM, nbroad, gamma_L, n_air, [], [], [], T, P)
    else:
        gamma = np.array([])

    if profile_label == 'Doppler':
        alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
        return cross_section_Doppler(wn_grid, v, alpha, coef, cutoff)
    elif profile_label == 'Gaussian':
        alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
        return cross_section_Doppler(wn_grid, v, alpha, coef, cutoff)
    elif profile_label == 'Lorentzian':
        return cross_section_Lorentzian(wn_grid, v, gamma, coef, cutoff)
    elif profile_label == 'SciPy Voigt':
        alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
        sigma = Gaussian_standard_deviation(alpha)
        return cross_section_SciPyVoigt(wn_grid, v, sigma, gamma, coef, cutoff)
    elif profile_label == 'SciPy wofz Voigt':
        alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
        sigma = Gaussian_standard_deviation(alpha)
        return cross_section_SciPyWofzVoigt(wn_grid, v, sigma, gamma, coef, cutoff)
    elif profile_label == 'Humlicek Voigt':
        alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
        return cross_section_HumlicekVoigt(wn_grid, v, alpha, gamma, coef, cutoff)
    elif profile_label == 'Thompson pseudo-Voigt':
        alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
        eta = PseudoThompsonVoigt(alpha, gamma)
        return cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, coef, cutoff)
    elif profile_label == 'Kielkopf pseudo-Voigt':
        alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
        eta = PseudoKielkopfVoigt(alpha, gamma)
        return cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, coef, cutoff)
    elif profile_label == 'Olivero pseudo-Voigt':
        alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
        eta = PseudoOliveroVoigt(alpha, gamma)
        return cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, coef, cutoff)
    elif profile_label == 'Liu-Lin pseudo-Voigt':
        alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
        eta = PseudoLiuLinVoigt(alpha, gamma)
        return cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, coef, cutoff)
    elif profile_label == 'Rocco pseudo-Voigt':
        alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
        eta = PseudoRoccoVoigt(alpha, gamma)
        return cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, coef, cutoff)
    elif profile_label == 'Binned Doppler':
        alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
        return cross_section_BinnedGaussian(wn_grid, v, alpha, coef, cutoff)
    elif profile_label == 'Binned Gaussion':
        alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
        return cross_section_BinnedGaussian(wn_grid, v, alpha, coef, cutoff)
    elif profile_label == 'Binned Lorentzian':
        return cross_section_BinnedLorentzian(wn_grid, v, gamma, coef, cutoff)
    elif profile_label == 'Binned Voigt':
        alpha = DopplerHWHM_alpha(num, alpha_HWHM, v, T)
        sigma = Gaussian_standard_deviation(alpha)
        return cross_section_BinnedVoigt(wn_grid, v, sigma, gamma, coef, cutoff)
    raise ValueError('Please choose line profile from the list.')


def save_exomolhr_cross_section(exomolhr_df, T_list, P_list, Tvib_list, Trot_list):
    """
    Calculate and save cross sections for ExoMolHR line lists.
    """
    from pyexocross.core import (
        NLTEMethod,
        abs_emi,
        cutoff,
        UncFilter,
        min_wnl,
        max_wnl,
        wn_grid,
        profile,
        read_path,
        database,
    )

    print('Calculating cross sections ...')
    t = Timer()
    t.start()

    print('Preparing ExoMolHR line list data ONCE for all temperatures...')
    _, _, _, _, _, Q_arr, Evibp, Erotp, Evibpp, Erotpp = process_exomolhr_linelist_Q(
        exomolhr_df,
        T_list,
        Tvib_list,
        Trot_list,
    )
    profile_label = line_profile(profile)
    pressure_dependent = profile_label not in ['Gaussian', 'Doppler', 'Binned Doppler', 'Binned Gaussion']

    n_temps = get_ntemp(NLTEMethod, T_list, Trot_list)
    if isinstance(P_list, (list, np.ndarray)) and len(P_list) > 0:
        first_item = P_list[0]
        if isinstance(first_item, (list, np.ndarray)):
            P_per_temp = [list(P_list[i]) if i < len(P_list) else list(P_list[0]) for i in range(n_temps)]
        else:
            P_values = list(P_list)
            P_per_temp = [P_values[:] for _ in range(n_temps)]
    else:
        P_single = P_list if isinstance(P_list, (int, float)) else (P_list[0] if len(P_list) > 0 else 1.0)
        P_per_temp = [[P_single] for _ in range(n_temps)]

    broad, ratio, _, _ = read_broad(read_path)
    print_xsec_info(profile_label, cutoff, UncFilter, min_wnl, max_wnl,
                    'cm⁻¹', 'cm⁻¹/(molecule cm⁻²)', broad, ratio)

    any_results = False
    xsec_file_count = 0
    for temp_idx in log_tqdm(range(n_temps), desc='\nProcessing cross sections'):
        T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
        pressures = P_per_temp[temp_idx] if pressure_dependent else [P_per_temp[temp_idx][0]]
        Q = Q_arr[temp_idx]
        for P in pressures:
            xsec = process_exomolhr_cross_section_chunk(
                exomolhr_df,
                T,
                P,
                Q,
                profile_label,
                Tvib,
                Trot,
                Evibp,
                Erotp,
                Evibpp,
                Erotpp,
            )
            if len(xsec) == 0 or np.all(xsec == 0):
                print(f'Warning: No cross sections found for T={T} K, P={P} bar. Skipping this combination.')
                continue

            any_results = True
            save_xsec_file_plot(wn_grid, xsec, database, profile_label, T, P, temp_idx, Tvib_list, Trot_list)
            xsec_file_count += 1
            print_T_Tvib_Trot_P_path_info(T, Tvib, Trot, P if pressure_dependent else None, abs_emi, NLTEMethod, 'Cross sections', None)

    print('\nTotal running time for cross sections:')
    t.end()
    print('\nFinished calculating cross sections!\n')

    if not any_results:
        raise ValueError("Empty result with the input filter values. Please type new filter values in the input file.")

    print(f'\nAll {xsec_file_count} cross sections files have been saved!\n')
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')
